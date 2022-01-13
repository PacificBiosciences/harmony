#include "HarmonySettings.hpp"
#include "LibraryInfo.hpp"
#include "SimpleBamParser.h"

#include <htslib/hts.h>
#include <pbbam/FastaReader.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbbam/PbbamVersion.h>
#include <pbcopper/cli2/CLI.h>
#include <pbcopper/cli2/internal/BuiltinOptions.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/MemoryConsumption.h>
#include <pbcopper/utility/PbcopperVersion.h>
#include <pbcopper/utility/Stopwatch.h>
#include <zlib.h>
#include <boost/version.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>

namespace PacBio {
namespace Harmony {

static constexpr std::array<char, 4> bases{'A', 'C', 'G', 'T'};

std::unordered_map<std::string, std::string> ReadRefs(const std::string& refFile)
{
    std::unordered_map<std::string, std::string> refs;
    PacBio::BAM::FastaReader fastaReader{refFile};
    PacBio::BAM::FastaSequence fasta;
    while (fastaReader.GetNext(fasta))
        refs.insert({fasta.Name(), fasta.Bases()});
    return refs;
}

void SetBamReaderDecompThreads(const int32_t numThreads)
{
    static constexpr char BAMREADER_ENV[] = "PB_BAMREADER_THREADS";
    const std::string decompThreads = std::to_string(numThreads);
    setenv(BAMREADER_ENV, decompThreads.c_str(), true);
}

std::string ParseAlignment(const BAM::BamRecord& record,
                           const std::unordered_map<std::string, std::string>& refs)
{
    std::ostringstream out;
    std::map<char, std::map<char, int32_t>> singleBase;
    std::map<char, int32_t> singleBaseDel;
    std::map<char, int32_t> allBaseDel;
    std::map<char, std::map<char, int32_t>> singleBaseIns;
    std::map<char, std::map<char, int32_t>> allBaseIns;
    int32_t ins = 0;
    int32_t del = 0;
    int32_t insEvents = 0;
    int32_t delEvents = 0;
    int32_t insMultiEvents = 0;
    int32_t delMultiEvents = 0;
    int32_t mismatch = 0;
    int32_t match = 0;
    const auto ref =
        refs.at(record.ReferenceName())
            .substr(record.ReferenceStart(), record.ReferenceEnd() - record.ReferenceStart());
    int32_t qryPos = 0;
    int32_t refPos = 0;
    const auto qry = record.Sequence(BAM::Orientation::GENOMIC);
    for (const auto& cigar : record.CigarData()) {
        int32_t len = cigar.Length();
        switch (cigar.Type()) {
            case BAM::CigarOperationType::INSERTION:
                ++singleBaseIns[ref[refPos]][qry[qryPos]];
                ++insEvents;
                if (len > 1) ++insMultiEvents;
                ins += len;
                qryPos += len;
                break;
            case BAM::CigarOperationType::DELETION:
                ++singleBaseDel[ref[refPos]];
                for (int32_t i = 0; i < len; ++i)
                    ++allBaseDel[ref[refPos + i]];
                ++delEvents;
                if (len > 1) ++delMultiEvents;
                del += len;
                refPos += len;
                break;
            case BAM::CigarOperationType::SEQUENCE_MISMATCH:
                for (int32_t i = 0; i < len; ++i) {
                    ++singleBase[ref[refPos + i]][qry[qryPos + i]];
                }
                mismatch += len;
                refPos += len;
                qryPos += len;
                break;
            case BAM::CigarOperationType::SEQUENCE_MATCH:
                for (int32_t i = 0; i < len; ++i) {
                    ++singleBase[ref[refPos + i]][qry[qryPos + i]];
                }
                refPos += len;
                qryPos += len;
                match += len;
                break;
            case BAM::CigarOperationType::SOFT_CLIP:
                qryPos += len;
                break;
            case BAM::CigarOperationType::ALIGNMENT_MATCH:
                PBLOG_FATAL << "UNSUPPORTED OPERATION: ALIGNMENT MATCH";
                std::exit(EXIT_FAILURE);
            case BAM::CigarOperationType::REFERENCE_SKIP:
                PBLOG_FATAL << "UNSUPPORTED OPERATION: REFERENCE SKIP";
                std::exit(EXIT_FAILURE);
            case BAM::CigarOperationType::HARD_CLIP:
                PBLOG_FATAL << "UNSUPPORTED OPERATION: HARD CLIP";
                std::exit(EXIT_FAILURE);
            case BAM::CigarOperationType::PADDING:
                PBLOG_FATAL << "UNSUPPORTED OPERATION: PADDING";
                std::exit(EXIT_FAILURE);
            case BAM::CigarOperationType::UNKNOWN_OP:
            default:
                PBLOG_FATAL << "UNKNOWN OP";
                std::exit(EXIT_FAILURE);
                break;
        }
    }
    const int32_t span = record.AlignedEnd() - record.AlignedStart();
    const int32_t nErr = ins + del + mismatch;
    const int32_t numAlignedBases = match + ins + mismatch;
    const double concordance = 1.0 - 1.0 * nErr / span;
    const int32_t qv = concordance == 1 ? 60 : (-10 * std::log10(1 - concordance));
    const int32_t numPasses = record.HasNumPasses() ? record.NumPasses() : -1;
    const int32_t ec = record.Impl().HasTag("ec") ? record.Impl().TagValue("ec").ToFloat() : -1;
    const std::string name = record.FullName();
    const double rq = record.ReadAccuracy();
    const int32_t seqlen = qry.size();

    out << name << ' ' << numPasses << ' ' << ec << ' ' << rq << ' ' << seqlen << ' '
        << numAlignedBases << ' ' << concordance << ' ' << qv << ' ' << match << ' ' << mismatch
        << ' ' << del << ' ' << ins << ' ' << delEvents << ' ' << insEvents << ' ' << delMultiEvents
        << ' ' << insMultiEvents;
    for (const auto& refBase : bases) {
        for (const auto& qryBase : bases) {
            if (singleBase.find(refBase) == singleBase.cend() ||
                singleBase.at(refBase).find(qryBase) == singleBase.at(refBase).cend()) {
                out << ' ' << 0;
            } else {
                out << ' ' << singleBase.at(refBase).at(qryBase);
            }
        }
    }
    for (const auto& refBase : bases) {
        for (const auto& qryBase : bases) {
            if (singleBaseIns.find(refBase) == singleBaseIns.cend() ||
                singleBaseIns.at(refBase).find(qryBase) == singleBaseIns.at(refBase).cend()) {
                out << ' ' << 0;
            } else {
                out << ' ' << singleBaseIns.at(refBase).at(qryBase);
            }
        }
    }
    for (const auto& refBase : bases) {
        if (singleBaseDel.find(refBase) == singleBaseDel.cend()) {
            out << ' ' << 0;
        } else {
            out << ' ' << singleBaseDel.at(refBase);
        }
    }
    for (const auto& refBase : bases) {
        for (const auto& qryBase : bases) {
            if (allBaseIns.find(refBase) == allBaseIns.cend() ||
                allBaseIns.at(refBase).find(qryBase) == allBaseIns.at(refBase).cend()) {
                out << ' ' << 0;
            } else {
                out << ' ' << allBaseIns.at(refBase).at(qryBase);
            }
        }
    }
    for (const auto& refBase : bases) {
        if (allBaseDel.find(refBase) == allBaseDel.cend()) {
            out << ' ' << 0;
        } else {
            out << ' ' << allBaseDel.at(refBase);
        }
    }
    out << '\n';
    return out.str();
}

int RunnerSubroutine(const CLI_v2::Results& options)
{
    Utility::Stopwatch globalTimer;
    HarmonySettings settings{options};
    SetBamReaderDecompThreads(settings.NumThreads);

    const std::string alnFile{settings.FileNames[0]};
    const std::string refFile{settings.FileNames[1]};

    std::unique_ptr<ReaderBase> alnReader = SimpleBamParser::BamQuery(alnFile, settings.Region);
    PBLOG_INFO << "Start reading reference";
    std::unordered_map<std::string, std::string> refs = ReadRefs(refFile);
    // std::unordered_map<std::string, std::vector<bool>> hps;
    // for (const auto& r : refs) {
    //     const auto& singleRef = r.second;
    //     const int32_t numBases = singleRef.size();
    //     std::vector<bool> hp(numBases, false);
    //     for (int32_t i = 0; i < numBases - 1; ++i) {
    //         if (singleRef[i] == singleRef[i + 1]) {
    //             hp[i] = true;
    //             hp[i + 1] = true;
    //         }
    //     }
    //     hps.insert({r.first, std::move(hp)});
    // }
    PBLOG_INFO << "Finished reading reference";

    std::array<char, 4> bases = {'A', 'C', 'G', 'T'};
    BAM::BamRecord record;

    std::ofstream outputFile{settings.FileNames[2]};

    outputFile << "name" << ' ' << "passes" << ' ' << "ec" << ' ' << "rq" << ' ' << "seqlen" << ' '
               << "alnlen" << ' ' << "concordance" << ' ' << "qv" << ' ' << "match" << ' '
               << "mismatch" << ' ' << "del" << ' ' << "ins" << ' ' << "del_events" << ' '
               << "ins_events" << ' ' << "del_multi_events" << ' ' << "ins_multi_events";
    for (const auto& refBase : bases) {
        for (const auto& qryBase : bases) {
            outputFile << " sub_" << refBase << qryBase;
        }
    }
    for (const auto& refBase : bases) {
        for (const auto& qryBase : bases) {
            outputFile << " ins_single_" << refBase << qryBase;
        }
    }
    for (const auto& refBase : bases) {
        outputFile << " del_single_" << refBase;
    }
    for (const auto& refBase : bases) {
        for (const auto& qryBase : bases) {
            outputFile << " ins_all_" << refBase << qryBase;
        }
    }
    for (const auto& refBase : bases) {
        outputFile << " del_all_" << refBase;
    }
    outputFile << '\n';

    int32_t counter = 0;
    while (alnReader->GetNext(record)) {
        if (++counter % 1000 == 0) {
            PBLOG_INFO << counter;
        }
        outputFile << ParseAlignment(record, refs);
    }

    globalTimer.Freeze();
    PBLOG_INFO << "Run Time : " << globalTimer.ElapsedTime();
    PBLOG_INFO << "CPU Time : "
               << Utility::Stopwatch::PrettyPrintNanoseconds(
                      static_cast<int64_t>(Utility::Stopwatch::CpuTime() * 1000 * 1000 * 1000));

    int64_t const peakRss = PacBio::Utility::MemoryConsumption::PeakRss();
    double const peakRssGb = peakRss / 1024.0 / 1024.0 / 1024.0;
    PBLOG_INFO << "Peak RSS : " << std::fixed << std::setprecision(3) << peakRssGb << " GB";

    return EXIT_SUCCESS;
}
}  // namespace Harmony
}  // namespace PacBio

int main(int argc, char* argv[])
{
    return PacBio::CLI_v2::Run(argc, argv, PacBio::Harmony::HarmonySettings::CreateCLI(),
                               &PacBio::Harmony::RunnerSubroutine);
}
