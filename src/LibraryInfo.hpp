#pragma once

#include <pbcopper/library/Bundle.h>
#include <pbcopper/library/Info.h>

namespace PacBio {
namespace Harmony {

///
/// \return Dummy library info (e.g. name, version)
///
Library::Info LibraryInfo();

///
/// \returns bundle of Dummy library info, plus its dependencies
///
Library::Bundle LibraryBundle();

}  // namespace Harmony
}  // namespace PacBio
