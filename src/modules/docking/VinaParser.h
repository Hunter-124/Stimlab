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

namespace stimlab::docking {

// Parse a Vina/smina output PDBQT string into ranked poses. `reference` supplies
// the per-atom z, the bond list and heavyCount that the coordinate-only PDBQT
// lacks; parsed coordinates are matched positionally to the reference atom order
// (heavy atoms + polar H, exactly as written by writeRigidPdbqt / obabel). Poses
// are returned sorted by ascending affinity (most negative = strongest first).
// Non-finite affinities/coordinates are skipped. Returns empty on no MODELs.
std::vector<DockPose> parseVinaPdbqt(const std::string& output,
                                     const chem::Conformer& reference);

}  // namespace stimlab::docking
