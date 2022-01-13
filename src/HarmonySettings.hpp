#pragma once

#include <pbcopper/cli2/CLI.h>

#include <cstdint>
#include <string>
#include <vector>

namespace PacBio {
namespace Harmony {

struct HarmonySettings
{
    const std::string CLI;
    const std::string LogFile;
    const std::vector<std::string> FileNames;
    const std::string Region;
    const int32_t NumThreads;
    const bool ExtendedMatrics;

    HarmonySettings(const PacBio::CLI_v2::Results& options);

    static CLI_v2::Interface CreateCLI();
};
}  // namespace Harmony
}  // namespace PacBio
