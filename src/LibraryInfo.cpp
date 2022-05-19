#include "LibraryInfo.hpp"

#include <LibraryGitHash.hpp>
#include <LibraryVersion.hpp>
#include <string>

namespace PacBio {
namespace Harmony {

Library::Bundle LibraryBundle()
{
    Library::Bundle bundle{LibraryInfo(), {}};
    return bundle;
}

Library::Info LibraryInfo()
{
    return {"Harmony", Harmony::ReleaseVersion, Harmony::LibraryGitSha1};
}

}  // namespace Harmony
}  // namespace PacBio
