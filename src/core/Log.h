// core/Log.h - logging bootstrap. Call init() once at startup.
// Sets up a colored console sink + a rotating file sink under AppPaths::logs(),
// and (on Windows) installs an unhandled-exception filter that writes a minidump.
#pragma once

#include <spdlog/spdlog.h>  // re-exported so call sites can use spdlog::info(...) etc.

namespace stimlab::log {

void init(bool installCrashHandler = true);
void shutdown();

}  // namespace stimlab::log
