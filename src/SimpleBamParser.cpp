// Author: Armin Töpfer

#include <pbcopper/utility/FileUtils.h>

#include "SimpleBamParser.h"

namespace PacBio {
BaiReader::BaiReader(const BAM::GenomicInterval& interval, const PacBio::BAM::DataSet& dataset)
    : query_{interval, dataset}
{}

bool BaiReader::GetNext(BAM::BamRecord& record) { return query_.GetNext(record); }

AlignedCollator::AlignedCollator(std::vector<std::unique_ptr<BAM::BamReader>> readers)
{
    for (auto&& reader : readers) {
        auto item = BAM::internal::CompositeMergeItem{std::move(reader)};
        if (item.reader->GetNext(item.record)) mergeItems_.push_back(std::move(item));
    }
    UpdateSort();
}

void AlignedCollator::UpdateSort()
{
    std::sort(mergeItems_.begin(), mergeItems_.end(), PacBio::BAM::PositionSorter{});
}

bool AlignedCollator::GetNext(BAM::BamRecord& record)
{
    // nothing left to read
    if (mergeItems_.empty()) return false;

    // non-destructive 'pop' of first item from queue
    auto firstIter = mergeItems_.begin();
    auto firstItem = BAM::internal::CompositeMergeItem{std::move(firstIter->reader),
                                                       std::move(firstIter->record)};
    mergeItems_.pop_front();

    // store its record in our output record
    std::swap(record, firstItem.record);

    // try fetch 'next' from first item's reader
    // if successful, re-insert it into container & re-sort on our new values
    // otherwise, this item will go out of scope & reader destroyed
    if (firstItem.reader->GetNext(firstItem.record)) {
        mergeItems_.push_front(std::move(firstItem));
        UpdateSort();
    }

    // return success
    return true;
}

std::vector<std::unique_ptr<BAM::BamReader>> SimpleBamParser::GetBamReaders(
    const std::string& filePath, const BAM::PbiFilter& filter)
{
    BAM::DataSet ds(filePath);
    std::vector<std::string> inputFilenames_;
    const auto& bamFiles = ds.BamFiles();
    inputFilenames_.reserve(bamFiles.size());
    for (const auto& file : bamFiles)
        inputFilenames_.push_back(file.Filename());

    if (inputFilenames_.empty())
        throw std::runtime_error("no input filenames provided to BamFileMerger");

    // attempt open input files
    std::vector<std::unique_ptr<BAM::BamReader>> readers;
    readers.reserve(inputFilenames_.size());
    for (const auto& fn : inputFilenames_) {
        if (filter.IsEmpty())
            readers.emplace_back(new BAM::BamReader(fn));
        else
            readers.emplace_back(new BAM::PbiIndexedBamReader(filter, fn));
    }
    return readers;
}

std::unique_ptr<ReaderBase> SimpleBamParser::BamQuery(const std::string& filePath,
                                                      const std::string& userFilters)
{
    if (!Utility::FileExists(filePath)) {
        PBLOG_FATAL << "Could not open input file " << filePath;
        std::exit(EXIT_FAILURE);
    }
    using namespace BAM;
    if (userFilters.empty()) return BamQuery(filePath);

    BAM::DataSet ds(filePath);
    int bamCounter = 0;
    int pbiCounter = 0;
    int baiCounter = 0;
    for (const auto& f : ds.BamFiles()) {
        ++bamCounter;
        pbiCounter += f.PacBioIndexExists();
        baiCounter += f.StandardIndexExists();
    }

    if (bamCounter == 0) {
        PBLOG_FATAL << "No input BAM files";
        std::exit(EXIT_FAILURE);
    }

    bool usePbi = bamCounter == pbiCounter;
    bool useBai = bamCounter == baiCounter;
    if (pbiCounter > 0 && baiCounter > 0) {
        PBLOG_WARN << "Both index files, pbi and bai are present.";
    }
    if (usePbi) {
        PBLOG_INFO << "Using PBI files for filtering";
        std::vector<PbiFilter> pbiFilters;
        std::vector<std::string> singleFilters;
        boost::split(singleFilters, userFilters, boost::is_any_of(";"));
        for (const auto& singleFilter : singleFilters) {
            std::vector<std::string> chr_pos;
            boost::split(chr_pos, singleFilter, boost::is_any_of(":"));

            PbiFilter refFilter = PbiReferenceNameFilter{chr_pos[0], Compare::EQUAL};
            if (chr_pos.size() == 1) {
                pbiFilters.emplace_back(std::move(refFilter));
            } else if (chr_pos.size() == 2) {
                std::vector<std::string> pos;
                boost::split(pos, chr_pos[1], boost::is_any_of("-"));
                if (pos.size() == 1) {
                    if (std::stoi(pos[0]) < 0) {
                        PBLOG_FATAL << "Reference position has to be non-negative.";
                        std::exit(EXIT_FAILURE);
                    }
                    boost::replace_all(pos[0], ",", "");
                    uint32_t singlePos = std::stoul(pos[0]);
                    pbiFilters.emplace_back(PbiFilter::Intersection(
                        {std::move(refFilter),
                         PbiReferenceEndFilter{singlePos, Compare::GREATER_THAN_EQUAL},
                         PbiReferenceStartFilter{singlePos, Compare::LESS_THAN_EQUAL}}));
                } else if (pos.size() == 2) {
                    boost::replace_all(pos[0], ",", "");
                    boost::replace_all(pos[1], ",", "");
                    if (std::stoi(pos[0]) < 0 || std::stoi(pos[1]) < 0) {
                        PBLOG_FATAL << "Reference position has to be non-negative.";
                        std::exit(EXIT_FAILURE);
                    }
                    uint32_t startPos = std::stoul(pos[0]);
                    uint32_t endPos = std::stoul(pos[1]);
                    pbiFilters.emplace_back(PbiFilter::Intersection(
                        {std::move(refFilter),
                         PbiReferenceEndFilter{startPos, Compare::GREATER_THAN_EQUAL},
                         PbiReferenceStartFilter{endPos, Compare::LESS_THAN_EQUAL}}));
                } else if (pos.size() > 2) {
                    PBLOG_FATAL << "Only two positions per filter allowed.";
                    std::exit(EXIT_FAILURE);
                }
            } else if (chr_pos.size() > 2) {
                PBLOG_FATAL << "Only one : per filter allowed.";
                std::exit(EXIT_FAILURE);
            }
        }

        const auto filter = BAM::PbiFilter::FromDataSet(ds);
        BAM::PbiFilter filters;
        if (filter.IsEmpty())
            filters = PbiFilter::Union(std::move(pbiFilters));
        else
            filters = PbiFilter::Intersection(
                {PbiFilter::Union(std::move(pbiFilters)), std::move(filter)});
        return std::make_unique<AlignedCollator>(GetBamReaders(filePath, filters));
    } else if (useBai) {
        PBLOG_INFO << "Using BAI files for filtering";
        std::string filterNoComma = userFilters;
        boost::replace_all(filterNoComma, ",", "");
        return std::make_unique<BaiReader>(filterNoComma, ds);
    }
    PBLOG_FATAL << "Number of index files does not match number of BAM files!";
    std::exit(EXIT_FAILURE);
}

std::unique_ptr<ReaderBase> SimpleBamParser::BamQuery(const std::string& filePath)
{
    BAM::DataSet ds(filePath);
    const auto filter = BAM::PbiFilter::FromDataSet(ds);
    return std::make_unique<AlignedCollator>(GetBamReaders(filePath, filter));
}

// # BamIO
BAM::BamHeader SimpleBamParser::ExtractHeader(const std::string& datasetPath)

{
    BAM::BamHeader header;
    const auto bamFiles = BAM::DataSet(datasetPath).BamFiles();
    if (bamFiles.empty()) {
        std::cerr << "No BAM files available for: " << datasetPath << std::endl;
        std::exit(EXIT_FAILURE);
    }

    header = bamFiles.front().Header();
    for (size_t i = 1; i < bamFiles.size(); ++i)
        header += bamFiles.at(i).Header();

    return header;
}
std::vector<BAM::ReadGroupInfo> SimpleBamParser::ExtractReadGroups(const std::string& datasetPath)
{
    std::vector<BAM::ReadGroupInfo> out;
    const auto bamFiles = BAM::DataSet(datasetPath).BamFiles();
    if (bamFiles.empty()) {
        std::cerr << "No BAM files available for: " << datasetPath << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < bamFiles.size(); ++i)
        for (const auto& rg : bamFiles.at(i).Header().ReadGroups())
            out.emplace_back(rg);

    return out;
}
}  // namespace PacBio
