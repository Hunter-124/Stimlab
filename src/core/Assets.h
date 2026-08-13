// core/Assets.h - locates the read-only asset tree shipped beside the executable.
//
// WHY this lives in core and not in packs: both the pack loader (compound and
// target catalogs) and the chem engine's rule packs (functional groups,
// structural alerts) have to find assets/, and biocad_chem cannot depend on
// biocad_packs - biocad_contracts already depends on biocad_chem, so that link
// would be a cycle. core is the one layer below both.
#pragma once

#include <filesystem>
#include <string_view>

namespace biocad::core {

// The shipped asset root: <exe dir>/assets when packaged, the source tree's
// assets/ when running from a dev build directory, else ./assets. Returns an
// empty path when no asset tree is found - callers must report that rather than
// silently behaving as if the data were empty.
std::filesystem::path assetRoot();

// assetRoot()/sub, or an empty path when assetRoot() is empty.
std::filesystem::path assetDir(std::string_view sub);

}  // namespace biocad::core
