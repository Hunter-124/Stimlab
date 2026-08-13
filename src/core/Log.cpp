#include "core/Log.h"

#include "core/AppPaths.h"

#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <memory>
#include <vector>

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#  include <dbghelp.h>
#  include <ctime>
#endif

namespace biocad::log {

namespace {

#if defined(_WIN32)
LONG WINAPI writeMinidump(EXCEPTION_POINTERS* info) {
    auto dir = AppPaths::instance().logs();
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);

    const std::time_t now = std::time(nullptr);
    wchar_t name[64];
    swprintf(name, 64, L"crash-%lld.dmp", static_cast<long long>(now));
    const auto path = dir / name;

    HANDLE file = CreateFileW(path.c_str(), GENERIC_WRITE, 0, nullptr,
                              CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (file != INVALID_HANDLE_VALUE) {
        MINIDUMP_EXCEPTION_INFORMATION mdei{};
        mdei.ThreadId = GetCurrentThreadId();
        mdei.ExceptionPointers = info;
        mdei.ClientPointers = FALSE;
        MiniDumpWriteDump(GetCurrentProcess(), GetCurrentProcessId(), file,
                          MiniDumpNormal, &mdei, nullptr, nullptr);
        CloseHandle(file);
    }
    spdlog::critical("Unhandled exception - wrote minidump");
    return EXCEPTION_EXECUTE_HANDLER;
}
#endif

}  // namespace

void init(bool installCrashHandler) {
    AppPaths::instance().ensureLayout();

    std::vector<spdlog::sink_ptr> sinks;
    sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());

    const auto logfile = (AppPaths::instance().logs() / "biocad.log").string();
    constexpr std::size_t kMaxSize = 5 * 1024 * 1024;  // 5 MiB
    constexpr std::size_t kMaxFiles = 3;
    sinks.push_back(std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
        logfile, kMaxSize, kMaxFiles));

    auto logger = std::make_shared<spdlog::logger>("biocad", sinks.begin(), sinks.end());
    logger->set_level(spdlog::level::debug);
    logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
    spdlog::set_default_logger(logger);
    spdlog::flush_on(spdlog::level::warn);

#if defined(_WIN32)
    if (installCrashHandler) {
        SetUnhandledExceptionFilter(writeMinidump);
    }
#else
    (void)installCrashHandler;
#endif

    spdlog::info("BioCAD logging initialized -> {}", logfile);
}

void shutdown() {
    spdlog::shutdown();
}

}  // namespace biocad::log
