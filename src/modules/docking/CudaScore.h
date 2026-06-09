// modules/docking/CudaScore.h - CUDA GPU pose scorer (CUDA-free interface).
//
// A thin, CUDA-free seam to the GPU kernel in CudaScore.cu: only this header is
// included by the C++ (MSVC) translation units, so the rest of the docking module never
// sees <cuda_runtime.h>. CudaScore.cu (compiled by nvcc) implements it. It scores many
// RIGID ligand poses against a rigid receptor with the AutoDock Vina inter-molecular
// scoring function (Trott & Olson 2010), entirely on the GPU.
//
// Atom typing is a 2-class simplification: 0 = hydrophobic (carbon/aromatic),
// 1 = polar (N/O, H-bond-capable). radii are Vina xs surface radii (Angstrom).
#pragma once

namespace stimlab::docking {

// Is a usable CUDA device present right now? false when this is a non-CUDA build (a
// stub is never linked - the whole file is compiled only with STIMLAB_ENABLE_CUDA) or
// when no NVIDIA GPU is visible at runtime.
bool cudaDockAvailable();

// Score nPoses rigid poses on the GPU. Receptor atoms are fixed; the ligand base
// coordinates (ligXYZ, nLig*3 floats [x0,y0,z0,...]) are rigidly transformed by each
// pose's 3x4 row-major transform (poseXform, nPoses*12 floats:
// [r00 r01 r02 tx  r10 r11 r12 ty  r20 r21 r22 tz]) on the GPU, then scored against the
// receptor. One Vina score per pose is written to outScores (more negative = stronger).
// Returns false on no CUDA device / allocation / launch error (caller falls back to CPU).
bool cudaScorePoses(const float* recXYZ, const int* recType, const float* recRadius, int nRec,
                    const float* ligXYZ, const int* ligType, const float* ligRadius, int nLig,
                    const float* poseXform, int nPoses, float* outScores);

}  // namespace stimlab::docking
