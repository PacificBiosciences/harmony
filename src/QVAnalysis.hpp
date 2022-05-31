#pragma once

#include <pbbam/BamRecord.h>

#include <atomic>
#include <utility>

namespace PacBio {
namespace Harmony {

class QVAnalysis
{

public:
    void ProcessRecord(const BAM::BamRecord& record);
    void ComputeEmpiricalQVs(const std::string &outputDestination) const;

private:
    std::array<std::pair<std::atomic_int64_t, std::atomic_int64_t>, 94> counts_ = {};
};
}  // namespace Harmony
}  // namespace PacBio
