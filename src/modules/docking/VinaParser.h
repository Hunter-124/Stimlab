// modules/docking/VinaParser.h - parse AutoDock Vina / smina output PDBQT.
//
// Vina writes a multi-MODEL PDBQT: each MODEL is one scored pose, prefixed by a
//   REMARK VINA RESULT:   <affinity_kcal/mol>  <rmsd_lb>  <rmsd_ub>
// line, followed by ATOM/HETATM records. This parser turns that text into ranked
// DockPose objects (most negative affinity first), each carrying a chem::Conformer
// rebuilt from the parsed coordinates. Because a Vina pose is just relocated
// coordinates of the SAME ligand we submitted, the element list (z), the bond list
// and heavyCount are taken from the reference conformer we embedded; only the xyz
// changes per pose. smina's output uses the identical REMARK convention, so the
// same parser serves both engines.
//
// SAFETY SCOPE: poses are docked binding geometries for a binding-affinity
// (pharmacology) readout. No synthesis/route content is parsed or produced.
#pragma once

#include <string>
#include <vector>

#include "chem/Embed3D.h"
#include "contracts/IDockingBackend.h"

namespace biocad::docking {

// Parse a Vina/smina output PDBQT string into ranked poses. `reference` supplies
// the per-atom z, the bond list and heavyCount that the coordinate-only PDBQT
// lacks. Poses are returned sorted by ascending affinity (strongest first).
//
// `serialToConf` maps each output ATOM (in file order) to its index in `reference`.
// The flexible writer reorders atoms into the torsion tree, so without this map an
// index-only match would scramble elements/bonds; pass the writer's serialToConf so
// docked coordinates are scattered back onto the original topology. When it is null
// or empty, coordinates are matched positionally (the legacy/rigid behaviour).
std::vector<DockPose> parseVinaPdbqt(const std::string& output,
                                     const chem::Conformer& reference,
                                     const std::vector<int>* serialToConf = nullptr);

}  // namespace biocad::docking
