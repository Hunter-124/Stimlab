// modules/docking/Provisioning.h - off-thread engine + receptor provisioning.
//
// A Provisioner runs the slow, best-effort work behind real docking - locating /
// downloading vina.exe and preparing the headline CNS receptor PDBQTs - on its own
// worker thread so the UI never blocks. The UI polls a thread-safe status snapshot
// each frame and offers a "Provision" button that flips allowDownload on. A
// locate-only probe (allowDownload=false) is cheap and touches no network, so it can
// run at startup to seed the initial status.
//
// SAFETY SCOPE: this provisions third-party docking BINARIES and prepares receptor
// structures for a binding-affinity (pharmacology) prediction. No synthesis content.
#pragma once

#include <atomic>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "contracts/IDockingBackend.h"
#include "core/Manifest.h"

namespace stimlab::docking {

// Record the currently-provisioned engine binaries + prepared receptors into
// %APPDATA%/StimLab/manifest.json (the self-heal source of truth). Called after a
// successful provision.
void writeRuntimeManifest();

// Self-heal on launch: load the manifest, delete any corrupt component so a later
// provision re-fetches it, and return the post-heal status. Safe with no manifest.
ManifestStatus selfHealManifest();

class Provisioner {
public:
    Provisioner() = default;
    ~Provisioner();
    Provisioner(const Provisioner&) = delete;
    Provisioner& operator=(const Provisioner&) = delete;

    // Spawn the worker once (no-op while one is already running). allowDownload=false
    // does a fast locate-only probe (no network); true fetches engine + receptors.
    // `headline` are the receptor targets to prepare (e.g. docking::headlinePresets()).
    void start(bool allowDownload, std::vector<ReceptorTarget> headline);

    [[nodiscard]] bool running() const { return running_.load(); }
    [[nodiscard]] bool everRun() const { return everRun_.load(); }
    [[nodiscard]] bool vinaReady() const { return vinaReady_.load(); }
    [[nodiscard]] bool obabelReady() const { return obabelReady_.load(); }
    [[nodiscard]] int  receptorsReady() const { return receptorsReady_.load(); }
    [[nodiscard]] int  receptorsTotal() const { return receptorsTotal_.load(); }
    [[nodiscard]] std::string status() const;  // thread-safe snapshot

private:
    void run(bool allowDownload, std::vector<ReceptorTarget> headline);
    void setStatus(std::string s);

    std::thread        worker_;
    std::atomic<bool>  running_{false};
    std::atomic<bool>  everRun_{false};
    std::atomic<bool>  vinaReady_{false};
    std::atomic<bool>  obabelReady_{false};
    std::atomic<int>   receptorsReady_{0};
    std::atomic<int>   receptorsTotal_{0};
    mutable std::mutex mu_;
    std::string        status_{"Docking engine not provisioned yet."};
};

}  // namespace stimlab::docking
