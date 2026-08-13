#include "modules/docking/Provisioning.h"

#include <filesystem>
#include <system_error>
#include <utility>

#include "core/AppPaths.h"
#include "core/Manifest.h"
#include "modules/docking/EngineLocator.h"
#include "modules/docking/ReceptorPrep.h"

namespace biocad::docking {

namespace fs = std::filesystem;

void writeRuntimeManifest() {
    Manifest m;
    if (auto vina = locateEngine(Engine::Vina)) m.addFromFile("engine:vina", "engine", *vina);
    if (auto obabel = locateEngine(Engine::Obabel)) m.addFromFile("engine:obabel", "engine", *obabel);
    std::error_code ec;
    const fs::path rdir = receptorsDir();
    if (fs::exists(rdir, ec) && fs::is_directory(rdir, ec)) {
        for (fs::directory_iterator it(rdir, ec), end; it != end && !ec; it.increment(ec)) {
            if (it->is_regular_file(ec) && it->path().extension() == ".pdbqt")
                m.addFromFile("receptor:" + it->path().stem().string(), "receptor", it->path());
        }
    }
    m.save(AppPaths::instance().manifest());
}

ManifestStatus selfHealManifest() {
    Manifest m = Manifest::load(AppPaths::instance().manifest());
    if (m.empty()) return {};
    m.heal();           // delete corrupt files so a later provision re-fetches them
    return m.verify();  // status after healing
}

Provisioner::~Provisioner() {
    if (worker_.joinable()) worker_.join();
}

std::string Provisioner::status() const {
    std::lock_guard<std::mutex> lk(mu_);
    return status_;
}

void Provisioner::setStatus(std::string s) {
    std::lock_guard<std::mutex> lk(mu_);
    status_ = std::move(s);
}

void Provisioner::start(bool allowDownload, std::vector<ReceptorTarget> headline) {
    bool expected = false;
    if (!running_.compare_exchange_strong(expected, true)) return;  // already running
    if (worker_.joinable()) worker_.join();  // reap a previously-finished worker
    everRun_.store(true);
    receptorsReady_.store(0);
    receptorsTotal_.store(static_cast<int>(headline.size()));
    worker_ = std::thread(&Provisioner::run, this, allowDownload, std::move(headline));
}

void Provisioner::run(bool allowDownload, std::vector<ReceptorTarget> headline) {
    // 1) Engine binaries.
    setStatus(allowDownload ? "Provisioning AutoDock Vina..." : "Checking for AutoDock Vina...");
    const auto vina = ensureVina(allowDownload);
    vinaReady_.store(vina.fetched);

    const auto obabel = ensureObabel();
    obabelReady_.store(obabel.fetched);

    // 2) Headline receptors.
    int ready = 0;
    for (size_t i = 0; i < headline.size(); ++i) {
        setStatus("Preparing receptor " + headline[i].id + " (" + std::to_string(i + 1) + "/" +
                  std::to_string(headline.size()) + ")...");
        const auto rec = ensureReceptor(headline[i], allowDownload);
        if (rec.ready) ++ready;
        receptorsReady_.store(ready);
    }

    // 3) Record what is now on disk into the manifest (self-heal source of truth).
    writeRuntimeManifest();

    // 4) Final summary.
    std::string s = vina.fetched ? "Vina ready. " : ("Vina unavailable: " + vina.note + " ");
    s += std::to_string(ready) + "/" + std::to_string(headline.size()) + " receptors prepared.";
    if (!allowDownload && (ready == 0 || !vina.fetched))
        s += " Click Provision to download the engine + receptors.";
    if (obabel.fetched) s += " (obabel: full protonation.)";
    setStatus(std::move(s));

    running_.store(false);
}

}  // namespace biocad::docking
