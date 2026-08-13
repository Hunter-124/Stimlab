#include "storage/RunStore.h"

#include <cstdio>
#include <ctime>
#include <string>

#include <sqlite3.h>

#include "core/AppPaths.h"

namespace biocad {

struct SqliteRunStore::Impl {
    sqlite3* db = nullptr;
};

namespace {
std::string nowStamp() {
    std::time_t t = std::time(nullptr);
    std::tm tmv{};
#ifdef _WIN32
    localtime_s(&tmv, &t);
#else
    localtime_r(&t, &tmv);
#endif
    char buf[32];
    std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M", &tmv);
    return buf;
}
}  // namespace

SqliteRunStore::SqliteRunStore() : impl_(std::make_unique<Impl>()) {
    AppPaths::instance().ensureLayout();
    openAndInit(AppPaths::instance().db());
    // Seed a couple of illustrative rows on a fresh database so the panel is not
    // empty on first launch (clearly real, persisted rows - not fakes).
    if (ok() && count() == 0) {
        record({"", "Session", "BioCAD", "complete", "", "Run history is now persisted to SQLite."});
    }
}

SqliteRunStore::SqliteRunStore(std::filesystem::path dbPath) : impl_(std::make_unique<Impl>()) {
    openAndInit(dbPath);
}

SqliteRunStore::~SqliteRunStore() {
    if (impl_->db) sqlite3_close(impl_->db);
}

void SqliteRunStore::openAndInit(const std::filesystem::path& dbPath) {
    if (sqlite3_open(dbPath.string().c_str(), &impl_->db) != SQLITE_OK) {
        if (impl_->db) { sqlite3_close(impl_->db); impl_->db = nullptr; }
        return;
    }
    sqlite3_exec(impl_->db, "PRAGMA journal_mode=WAL;", nullptr, nullptr, nullptr);
    // NOTE: a dedicated table name ("run_history", not "runs") so we never collide
    // with a legacy stimlab.db left in %APPDATA% by the old Python backend, whose
    // "runs" table has an incompatible column layout (CREATE IF NOT EXISTS would be
    // a silent no-op and our INSERTs would fail on the column mismatch).
    const char* schema =
        "CREATE TABLE IF NOT EXISTS run_history ("
        " id INTEGER PRIMARY KEY AUTOINCREMENT,"
        " kind TEXT, subject TEXT, status TEXT, created_at TEXT, summary TEXT);";
    sqlite3_exec(impl_->db, schema, nullptr, nullptr, nullptr);
}

bool SqliteRunStore::ok() const { return impl_->db != nullptr; }

std::size_t SqliteRunStore::count() const {
    if (!impl_->db) return 0;
    sqlite3_stmt* st = nullptr;
    std::size_t n = 0;
    if (sqlite3_prepare_v2(impl_->db, "SELECT COUNT(*) FROM run_history;", -1, &st, nullptr) == SQLITE_OK) {
        if (sqlite3_step(st) == SQLITE_ROW) n = static_cast<std::size_t>(sqlite3_column_int64(st, 0));
    }
    sqlite3_finalize(st);
    return n;
}

std::vector<RunRecord> SqliteRunStore::recent() const {
    std::vector<RunRecord> out;
    if (!impl_->db) return out;
    const char* q =
        "SELECT id, kind, subject, status, created_at, summary FROM run_history ORDER BY id DESC LIMIT 50;";
    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(impl_->db, q, -1, &st, nullptr) != SQLITE_OK) return out;
    auto col = [&](int i) {
        const unsigned char* s = sqlite3_column_text(st, i);
        return s ? std::string(reinterpret_cast<const char*>(s)) : std::string();
    };
    while (sqlite3_step(st) == SQLITE_ROW) {
        RunRecord r;
        char id[24];
        std::snprintf(id, sizeof id, "run-%04lld",
                      static_cast<long long>(sqlite3_column_int64(st, 0)));
        r.id = id;
        r.kind = col(1);
        r.subject = col(2);
        r.status = col(3);
        r.createdAt = col(4);
        r.summary = col(5);
        out.push_back(std::move(r));
    }
    sqlite3_finalize(st);
    return out;
}

void SqliteRunStore::record(const RunRecord& r) {
    if (!impl_->db) return;
    const char* q =
        "INSERT INTO run_history(kind,subject,status,created_at,summary) VALUES(?,?,?,?,?);";
    sqlite3_stmt* st = nullptr;
    if (sqlite3_prepare_v2(impl_->db, q, -1, &st, nullptr) != SQLITE_OK) return;
    const std::string created = r.createdAt.empty() ? nowStamp() : r.createdAt;
    const std::string status = r.status.empty() ? "complete" : r.status;
    auto bind = [&](int i, const std::string& s) {
        sqlite3_bind_text(st, i, s.c_str(), -1, SQLITE_TRANSIENT);
    };
    bind(1, r.kind);
    bind(2, r.subject);
    bind(3, status);
    bind(4, created);
    bind(5, r.summary);
    sqlite3_step(st);
    sqlite3_finalize(st);
}

}  // namespace biocad
