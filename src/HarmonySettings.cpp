#include "HarmonySettings.hpp"

#include <htslib/hts.h>
#include <pbbam/PbbamVersion.h>
#include <pbcopper/cli2/internal/BuiltinOptions.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/PbcopperVersion.h>
#include <zlib.h>

#include <boost/algorithm/string/replace.hpp>
#include <boost/version.hpp>
#include <iostream>

#include "LibraryInfo.hpp"

namespace PacBio {
namespace Harmony {
namespace OptionNames {
// clang-format off
const CLI_v2::Option Region {
R"({
    "names" : ["region"],
    "description" : "Genomic region",
    "type" : "string",
    "default" : ""
})"
};
const CLI_v2::Option ExtendedMatrics {
R"({
    "names" : ["e", "extended-metrics"],
    "description" : "Output extended metrics, not required for harmony plots",
    "type" : "bool"
})"
};
// clang-format on
}  // namespace OptionNames

HarmonySettings::HarmonySettings(const PacBio::CLI_v2::Results& options)
    : CLI(options.InputCommandLine())
    , LogFile(options[CLI_v2::Builtin::LogFile])
    , FileNames(options.PositionalArguments())
    , Region(options[OptionNames::Region])
    , NumThreads(options.NumThreads())
    , ExtendedMatrics(options[OptionNames::ExtendedMatrics])
{
    if (FileNames.size() > 4 || FileNames.size() < 2) {
        PBLOG_FATAL
            << "Please specify input alignment BAM file, optional reference FASTA file, and output "
               "harmony TSV file. Please see --help for more information.";
        std::exit(EXIT_FAILURE);
    }

    if ((ExtendedMatrics || !Region.empty()) && FileNames.size() != 3) {
        PBLOG_FATAL << "Please specify input alignment BAM file, reference FASTA file, and output "
                       "harmony TSV file. Please see --help for more information.";
        std::exit(EXIT_FAILURE);
    }
}

CLI_v2::Interface HarmonySettings::CreateCLI()
{
    static const std::string description{"Compute error profiles from alignments."};
    CLI_v2::Interface i{"harmony", description, Harmony::LibraryInfo().Release};

    Logging::LogConfig logConfig;
    logConfig.Header = "| ";
    logConfig.Delimiter = " | ";
    logConfig.Fields = Logging::LogField::TIMESTAMP | Logging::LogField::LOG_LEVEL;
    i.LogConfig(logConfig);

    const CLI_v2::PositionalArgument inputAlignFile{
        R"({
        "name" : "IN.aligned.bam",
        "description" : "Aligned BAM.",
        "type" : "file",
        "required" : false
    })"};
    const CLI_v2::PositionalArgument inputRefFile{
        R"({
        "name" : "IN.ref.fasta",
        "description" : "Reference FASTA.",
        "type" : "file",
        "required" : true
    })"};
    const CLI_v2::PositionalArgument outputHarmonyFile{
        R"({
        "name" : "OUT.harmony.txt",
        "description" : "Harmony TXT.",
        "type" : "file",
        "required" : true
    })"};
    i.AddPositionalArguments({inputAlignFile, inputRefFile, outputHarmonyFile});
    i.AddOption(OptionNames::Region);
    i.AddOption(OptionNames::ExtendedMatrics);

    const auto printVersion = [](const CLI_v2::Interface& interface) {
        const std::string harmonyVersion = []() {
            return Harmony::LibraryInfo().Release + " (commit " + Harmony::LibraryInfo().GitSha1 +
                   ')';
        }();
        const std::string pbbamVersion = []() { return BAM::LibraryFormattedVersion(); }();
        const std::string pbcopperVersion = []() {
            return Utility::LibraryVersionString() + " (commit " + Utility::LibraryGitSha1String() +
                   ')';
        }();
        const std::string boostVersion = []() {
            std::string v = BOOST_LIB_VERSION;
            boost::replace_all(v, "_", ".");
            return v;
        }();
        const std::string htslibVersion = []() { return std::string{hts_version()}; }();
        const std::string zlibVersion = []() { return std::string{ZLIB_VERSION}; }();

        std::cout << interface.ApplicationName() << " " << interface.ApplicationVersion() << '\n';
        std::cout << '\n';
        std::cout << "Using:\n";
        std::cout << "  harmony  : " << harmonyVersion << '\n';
        std::cout << "  pbbam    : " << pbbamVersion << '\n';
        std::cout << "  pbcopper : " << pbcopperVersion << '\n';
        std::cout << "  boost    : " << boostVersion << '\n';
        std::cout << "  htslib   : " << htslibVersion << '\n';
        std::cout << "  zlib     : " << zlibVersion << '\n';
    };
    i.RegisterVersionPrinter(printVersion);

    return i;
}
}  // namespace Harmony
}  // namespace PacBio
