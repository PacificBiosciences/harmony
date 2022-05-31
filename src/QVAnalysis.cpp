#include "QVAnalysis.hpp"

#include <pbcopper/logging/Logging.h>

#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <ostream>
#include <regex>
#include <stdexcept>

#include "../third-party/incbin/incbin.h"

namespace PacBio {
namespace Harmony {

INCTXT(QV_ANALYSIS, "qv-pred-emp.html");  //NOLINT(hicpp-no-assembler)

void QVAnalysis::ProcessRecord(const BAM::BamRecord& record)
{
    const auto& qvs = record.Qualities(BAM::Orientation::NATIVE);
    const auto& cigar = record.CigarData();

    int pos = 0;
    const auto getQV = [&]() {
        assert(pos < std::ssize(qvs));
        const int qv = qvs[pos++];
        assert(qv >= 0 && qv <= 93);
        return qv;
    };

    for (const auto cigarOp : cigar) {
        const int len = cigarOp.Length();
        const auto opType = cigarOp.Type();
        switch (opType) {
            case BAM::CigarOperationType::INSERTION:
            case BAM::CigarOperationType::SEQUENCE_MISMATCH:
                for (int i = 0; i < len; ++i) {
                    const int qv = getQV();
                    counts_[qv].second += 1;
                }
                break;
            case BAM::CigarOperationType::DELETION:
                break;
            case BAM::CigarOperationType::SEQUENCE_MATCH:
                for (int i = 0; i < len; ++i) {
                    const int qv = getQV();
                    counts_[qv].first += 1;
                }
                break;
            case BAM::CigarOperationType::SOFT_CLIP:
                pos += len;
                break;
            default:
                PBLOG_CRITICAL << record.FullName()
                               << ": Unsupported CIGAR operation: " << static_cast<int>(opType)
                               << " at pos=" << pos;
                throw std::runtime_error("Unsupported CIGAR Operation");
        }
    }
    assert(pos == ssize(qvs));
}

void QVAnalysis::ComputeEmpiricalQVs(const std::string& outputDestination) const
{
    std::stringstream ss;
    ss << "#PredictedQV,EmpiricalQV,BaseCount" << std::endl;
    int qvPredicted = 0;
    for (const auto& [hit, miss] : counts_) {
        const auto baseCount = hit + miss;
        double errorProb = baseCount > 0.0 ? static_cast<double>(miss) / baseCount : 0.0;
        int qvEmpirical =
            round(-10.0 * log10(std::max(std::numeric_limits<double>::epsilon(), errorProb)));
        ss << qvPredicted << "," << qvEmpirical << "," << baseCount << std::endl;
        ++qvPredicted;
    }

    if (outputDestination.ends_with(".html")) {
        std::ofstream fs(outputDestination.c_str());
        std::string htmlTemplate(gQV_ANALYSISData);
        fs << std::regex_replace(htmlTemplate, std::regex("const DATA = sample; \\/\\/ REPLACE"),
                                 "const DATA = `" + ss.str() + "`;");
    } else if (outputDestination == "-") {
        std::cout << ss.str();
    } else {
        // assume comma separated value file
        std::ofstream fs(outputDestination.c_str());
        fs << ss.str();
    }
}

}  // namespace Harmony
}  // namespace PacBio
