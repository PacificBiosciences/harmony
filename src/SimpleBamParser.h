// Author: Armin Töpfer
#pragma once

#include <pbbam/BamReader.h>
#include <pbbam/CompositeBamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/GenomicIntervalQuery.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiIndexedBamReader.h>
#include <pbbam/internal/QueryBase.h>
#include <pbcopper/logging/Logging.h>

#include <boost/algorithm/string.hpp>
#include <deque>
#include <memory>
#include <string>
#include <vector>

#include "pbbam/GenomicInterval.h"

namespace PacBio {

class ReaderBase
{
public:
    virtual ~ReaderBase(){};

    virtual bool GetNext(BAM::BamRecord&) = 0;
};

class BaiReader : public ReaderBase
{
public:
    BaiReader(const BAM::GenomicInterval& interval, const PacBio::BAM::DataSet& dataset);
    ~BaiReader() override = default;

    bool GetNext(BAM::BamRecord& record) override;

private:
    BAM::GenomicIntervalCompositeBamReader query_;
};

class AlignedCollator : public ReaderBase
{
public:
    explicit AlignedCollator(std::vector<std::unique_ptr<BAM::BamReader>> readers);
    ~AlignedCollator() override = default;

    bool GetNext(BAM::BamRecord& record) override;

private:
    void UpdateSort();

    std::deque<BAM::internal::CompositeMergeItem> mergeItems_;
};

struct SimpleBamParser
{
    static std::vector<std::unique_ptr<BAM::BamReader>> GetBamReaders(const std::string& filePath,
                                                                      const BAM::PbiFilter& filter);

    static std::unique_ptr<ReaderBase> BamQuery(const std::string& filePath,
                                                const std::string& userFilters);

    static std::unique_ptr<ReaderBase> BamQuery(const std::string& filePath);

    // # BamIO
    static BAM::BamHeader ExtractHeader(const std::string& datasetPath);
    static std::vector<BAM::ReadGroupInfo> ExtractReadGroups(const std::string& datasetPath);
};
}  // namespace PacBio
