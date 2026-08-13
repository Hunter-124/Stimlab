# Rendering: bonds the file does not contain, and the cartoon built from them

Phase 6 added two things the viewport did not have: a **bond graph for protein coordinates** and
a **generic indexed-mesh draw path** so a ribbon can be drawn at all. Before this, `MolViewport`
could only instance spheres and cylinders, and `bio::Structure` had no bonds, because a PDB or
mmCIF file does not contain any.

| Path | Contents |
| --- | --- |
| `assets/packs/structure/residue-templates.json` | Heavy-atom connectivity of the standard residues, by atom name |
| `src/bio/Connectivity.{h,cpp}` | Template application, polymer links, chain gaps, disulfides, diagnostics |
| `src/bio/Cartoon.{h,cpp}` | Guide points, spline frames, cross-sections, indexed mesh |
| `src/render/MolViewport.{h,cpp}` | `IndexedMesh` + `MeshVertexRgba` and the DX11 mesh pass |
| `src/ui/Panels.cpp` (`proteinStructure`) | The `Structure3D` cartoon / spheres-and-sticks toggle |

## Why a connectivity pack exists

`ATOM` records are positions and names. `CONECT` is optional, is only written for non-standard
groups, and is absent from most entries; `_struct_conn` describes links, not intra-residue
chemistry. So bonds have to come from the one piece of chemical identity the file does state -
the residue name and the atom names - and that is what the pack is: 31 residue templates (20
amino acids, the HIS tautomer and CHARMM/AMBER spellings as aliases, MSE, SEC, and A/C/G/U plus
DA/DC/DG/DT/DU), 395 named bond pairs with chemical bond orders.

Three rules make the result auditable rather than plausible:

- **An unrecognised residue is reported, never guessed.** `ConnectivityDiagnostics::
  unknownResidues` carries the name, occurrence count and atom count, and a warning string names
  it. A distance fallback exists (`ConnectivityOptions::inferByDistance`) but is **off by
  default** and marks every bond it creates `BondKind::DistanceInferred`.
- **A chain break is a gap, not a bond.** `C(i)-N(i+1)` is only created inside a 1.10-1.60 A
  window (`O3'-P` inside 1.30-1.90 A). Outside it, a `ChainGap` records the measured distance
  and no bond is made. A 12 A "peptide bond" across a disordered loop is a fabricated covalent
  link between residues that are nowhere near each other.
- **Alternate locations do not cross-bond.** Two atoms bond only when their altLoc indicators are
  compatible (blank, or the same letter), so the A and B copies of a disordered side chain are
  not fused into a branch that exists in no conformer.

Disulfides are found from `SG..SG` geometry (1.60-2.50 A) because `SSBOND` is optional, and are
reported as their own `BondKind` so they are never confused with a template bond.

## The cartoon pipeline

1. **Guide points.** One per residue with a CA. The ribbon normal is the carbonyl direction
   `O - C` (the peptide plane); when the carbonyl is absent the normal falls back to the plane of
   three consecutive CAs and the guide is flagged `normalFromGeometry`. Helix guides are moved
   `standardShift` of the way towards the midpoint of their neighbours, because a helical CA sits
   off the axis and the raw CA path is itself a coil.
2. **Spline.** A cardinal/Catmull-Rom cubic Hermite through the guide points. `h00(0) = 1` and
   every other basis function vanishes at `u = 0`, so a knot sample **is** its guide point - the
   test measures the residual as exactly `0.000e+00 A`. Phantom control points one residue
   spacing past each terminus (`overhangFactor`) give the end segments a defined tangent.
3. **Frames.** Tangent from the analytic derivative; normal interpolated between the two guide
   normals and then Gram-Schmidt orthogonalised against the tangent; binormal from the cross
   product. Worst orthonormality residual measured over 193 frames: `4.996e-16`.
4. **Cross-sections.** Ellipse for helix, flat four-sided rectangle for strand with an arrowhead
   over the strand's last residue, small circle for coil, plus a flat cap at each chain end.
5. **Extrusion.** Consecutive rings stitched into a triangle list with per-vertex positions,
   normals and a palette colour index.

### The constants

The values are Mol\*'s cartoon defaults. What each one does *here*:

| Constant | Value | Role in this implementation |
| --- | --- | --- |
| standard tension | 0.5 | Hermite tangent scale for coil and strand; 0.5 is classic uniform Catmull-Rom |
| helix tension | 0.9 | The same scale for helices: longer tangents keep the ribbon on the helical path instead of cutting each corner |
| standard shift | 0.5 | Fraction of the way a helix guide point moves towards the midpoint of its neighbours |
| overhang factor | 2 | Phantom control point length past a terminus, in half-residue spacings (2 = one full spacing) |
| size factor | 0.2 A | Base half-thickness of every cross-section |
| aspect ratio | 5 | Width / thickness of a flat cross-section, so a strand ribbon is 2.0 A wide and 0.4 A thick |
| arrow factor | 1.5 | Width multiplier at the base of a strand arrowhead, tapering to 0 at the tip |
| linear segments | 8 | Spline samples per residue |
| radial segments | 16 | Points around one cross-section ring (rounded down to a multiple of 4 so a rectangle and an ellipse emit the same count) |

### The 180-degree flip check, and why it exists

Successive peptide planes along a chain point in **alternating** directions: in an extended
strand the carbonyl of residue *i* points one way and residue *i+1*'s points nearly the opposite
way. Interpolating those raw normals rotates the ribbon by roughly half a turn per residue, and a
flat sheet renders as a corkscrew. `guidePoints()` therefore compares each normal with its
already-corrected predecessor and negates it when the dot product is negative.

This is not asserted, it is **measured**. `bio::totalTwistDegrees()` sums the absolute rotation of
the ribbon normal about the local tangent over consecutive samples, and the test runs the same
chain with `CartoonOptions::flipCheck` on and off. On the 10-residue ideal strand of
`tests/fixtures/ideal_helix_strand.pdb`:

| flipCheck | total twist | per residue |
| --- | --- | --- |
| off | 1660.71 deg | **184.5 deg** |
| on | 56.89 deg | **6.3 deg** |

184.5 deg/residue is the predicted corkscrew: half a turn per residue, exactly the alternation of
the peptide planes. Over the whole 25-residue chain the numbers are 2839.74 deg (118.3/residue)
off versus 1008.24 deg (42.0/residue) on; the helix contributes genuine twist even with the check
enabled, because its carbonyls really do rotate about the helix axis, which is why the strong
per-residue claim is made on the strand.

### Measured triangle counts

8 samples per residue x 16 radial segments x 2 triangles = **256 triangles per residue interval**,
plus a 16-triangle cap at each end. Measured on the fixture:

| Region | Residues | Vertices | Triangles | Triangles / residue |
| --- | --- | --- | --- | --- |
| ideal alpha helix | 12 | 1426 | 2848 | 237.3 |
| ideal beta strand | 10 | 1170 | 2336 | 233.6 |
| whole chain (helix + coil + strand) | 25 | 3090 | 6176 | 247.0 |

The per-residue figure approaches 256 from below because the count is per residue *interval*
(`n-1` of them). Extrapolated: a 300-residue chain is **76,576 triangles**, which is the ~77k
figure the plan quotes and is trivial for DX11.

## The DX11 mesh path, and what was NOT executed

`MolScene` gained an optional `IndexedMesh` drawn after the sphere and cylinder passes:
`MeshVertexRgba` (position, normal, packed `0xAABBGGRR` colour; 28 bytes, offsets 0/12/24),
DYNAMIC vertex and index buffers that grow by doubling exactly like the existing instance
buffers, an input layout of `POSITION` / `NORMAL` / `COLOR` (`R8G8B8A8_UNORM`), a `kMeshShader`
compiled at runtime beside the atom and bond shaders, and one `DrawIndexed`.

The ribbon pass uses a second rasterizer state with `CullMode = D3D11_CULL_NONE` and a pixel
shader that flips the normal towards the viewer. That is a deliberate choice, not laziness: a
cross-section's triangle winding depends on the sign of the local frame, which the geometry stage
does not normalise, and a flat ribbon is legitimately viewed from behind. Culling it would drop
faces on half the chain. The previous state is restored immediately after the draw so the next
frame's sphere and cylinder passes still get back-face culling.

**Honest statement of verification.** This repository targets Windows and `src/render/
MolViewport.cpp` includes `<windows.h>`, `<d3d11.h>` and `<d3dcompiler.h>`; it **cannot be
compiled, let alone executed, on the Linux workstation this phase was developed on.** No pixel
produced by the mesh pass has been seen. What was actually verified:

- **The geometry, by execution.** `tests/test_bio_cartoon.cpp` and
  `tests/test_bio_connectivity.cpp` were built against a real amalgamated Catch2 and run: 13 test
  cases, 38041 assertions, all passing. Every number in the tables above is that run's output.
- **The vertex layout, by execution.** `sizeof(MeshVertexRgba) == 28` with `nx` at 12 and `rgba`
  at 24 is `static_assert`ed and was compiled and run, and those are exactly the byte offsets
  declared in the `meshElems` input layout. A mismatch here is the classic silent garbage-geometry
  bug, and it is the part of the GPU contract that can be checked off-GPU.
- **`computeSceneBounds()`, by execution.** Extracted verbatim from `MolViewport.cpp` into a
  throwaway translation unit and run: a mesh-only scene and a mesh-plus-atom scene both frame
  correctly, so the camera fits a ribbon that contains no `AtomInst` at all.
- **`MolViewport.h` and `src/ui/Panels.cpp`, by compilation.** Both were syntax-checked with the
  real Dear ImGui headers; the panel's toggle, scene caching and diagnostics compile.
- **The DX11 calls themselves, by line-by-line review** against the surrounding code: buffer
  descriptors copy the existing `ensureAtomInstCapacity` pattern (`DYNAMIC` + `CPU_ACCESS_WRITE`,
  doubling, capacity zeroed on failure), the map/unmap pairs match the existing
  `MAP_WRITE_DISCARD` usage, the new shader reads the same `cbuffer CB` register `b0` layout, and
  every new device object is released in `releaseAll()` with the same `safeRelease` helper.

What remains unverified on this machine, and must be checked on the first Windows run: that
`D3DCompile` accepts `kMeshShader`, that `CreateInputLayout` validates against its vertex
shader signature, and that the ribbon appears with helices and strands distinguishable and no
corkscrewing - the last being the visual counterpart of the twist number above. The capture route
is `BioCAD.exe --shot docs/media/structure3d.png --shot-panel Structure3D`.
