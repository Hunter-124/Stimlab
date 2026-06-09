// modules/docking/CudaScore.cu - GPU pose scorer (compiled by nvcc).
//
// Implements CudaScore.h: scores many rigid ligand poses against a rigid receptor on
// the GPU using the AutoDock Vina inter-molecular scoring function (Trott & Olson,
// J Comput Chem 2010). One CUDA thread scores one pose: it rigidly transforms each
// ligand atom by the pose's 3x4 matrix and sums the Vina pair terms over receptor atoms
// within an 8 A cutoff. The weights/terms below are Vina's published intermolecular
// function; the 2-class atom typing (hydrophobic vs polar) is a documented
// simplification (see CudaBackend.cpp). This is a rigid-body search - ligand torsions
// are NOT varied here (that is AutoDock Vina's flexible CPU job).
#include "modules/docking/CudaScore.h"

#include <cuda_runtime.h>
#include <math.h>

namespace stimlab::docking {
namespace {

// Vina intermolecular pair energy for a surface distance `surf` = r - (R_i + R_j).
// t* are the 2-class atom types: 0 = hydrophobic, 1 = polar (H-bond-capable).
__device__ inline float vinaPair(float surf, int t1, int t2) {
    const float g1  = __expf(-(surf * surf) * 4.0f);          // gauss1: exp(-(surf/0.5)^2)
    const float d3  = surf - 3.0f;
    const float g2  = __expf(-(d3 * d3) * 0.125f);            // gauss2: exp(-((surf-3)/2)^2)
    const float rep = surf < 0.0f ? surf * surf : 0.0f;       // repulsion (only when overlapping)
    float hyd = 0.0f;                                          // hydrophobic ramp 1@0.5 -> 0@1.5
    if (t1 == 0 && t2 == 0) {
        if (surf < 0.5f) hyd = 1.0f;
        else if (surf < 1.5f) hyd = 1.5f - surf;
    }
    float hb = 0.0f;                                           // H-bond ramp 1@-0.7 -> 0@0
    if (t1 == 1 && t2 == 1) {
        if (surf < -0.7f) hb = 1.0f;
        else if (surf < 0.0f) hb = -surf * (1.0f / 0.7f);
    }
    return -0.035579f * g1 - 0.005156f * g2 + 0.840245f * rep - 0.035069f * hyd - 0.587439f * hb;
}

__global__ void scoreKernel(const float* rec, const int* recType, const float* recRad, int nRec,
                            const float* lig, const int* ligType, const float* ligRad, int nLig,
                            const float* xform, int nPoses, float* out) {
    const int p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= nPoses) return;
    const float* M = xform + p * 12;  // 3x4 row-major transform for this pose
    float score = 0.0f;
    for (int i = 0; i < nLig; ++i) {
        const float lx = lig[i * 3 + 0], ly = lig[i * 3 + 1], lz = lig[i * 3 + 2];
        const float x = M[0] * lx + M[1] * ly + M[2] * lz + M[3];
        const float y = M[4] * lx + M[5] * ly + M[6] * lz + M[7];
        const float z = M[8] * lx + M[9] * ly + M[10] * lz + M[11];
        const int   ti = ligType[i];
        const float ri = ligRad[i];
        for (int j = 0; j < nRec; ++j) {
            const float dx = x - rec[j * 3 + 0];
            const float dy = y - rec[j * 3 + 1];
            const float dz = z - rec[j * 3 + 2];
            const float r = sqrtf(dx * dx + dy * dy + dz * dz);
            if (r >= 8.0f) continue;  // pair cutoff
            score += vinaPair(r - (ri + recRad[j]), ti, recType[j]);
        }
    }
    out[p] = score;
}

template <class T>
T* devAlloc(size_t n) {
    T* p = nullptr;
    if (cudaMalloc(reinterpret_cast<void**>(&p), n * sizeof(T)) != cudaSuccess) return nullptr;
    return p;
}

}  // namespace

bool cudaDockAvailable() {
    int n = 0;
    return cudaGetDeviceCount(&n) == cudaSuccess && n > 0;
}

bool cudaScorePoses(const float* recXYZ, const int* recType, const float* recRadius, int nRec,
                    const float* ligXYZ, const int* ligType, const float* ligRadius, int nLig,
                    const float* poseXform, int nPoses, float* outScores) {
    int devs = 0;
    if (cudaGetDeviceCount(&devs) != cudaSuccess || devs == 0) return false;
    if (nRec <= 0 || nLig <= 0 || nPoses <= 0) return false;

    float* dRec  = devAlloc<float>(static_cast<size_t>(nRec) * 3);
    int*   dRecT = devAlloc<int>(nRec);
    float* dRecR = devAlloc<float>(nRec);
    float* dLig  = devAlloc<float>(static_cast<size_t>(nLig) * 3);
    int*   dLigT = devAlloc<int>(nLig);
    float* dLigR = devAlloc<float>(nLig);
    float* dX    = devAlloc<float>(static_cast<size_t>(nPoses) * 12);
    float* dOut  = devAlloc<float>(nPoses);

    bool ok = dRec && dRecT && dRecR && dLig && dLigT && dLigR && dX && dOut;
    if (ok) {
        cudaMemcpy(dRec, recXYZ, static_cast<size_t>(nRec) * 3 * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(dRecT, recType, nRec * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dRecR, recRadius, nRec * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(dLig, ligXYZ, static_cast<size_t>(nLig) * 3 * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(dLigT, ligType, nLig * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dLigR, ligRadius, nLig * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(dX, poseXform, static_cast<size_t>(nPoses) * 12 * sizeof(float), cudaMemcpyHostToDevice);

        const int threads = 128;
        const int blocks = (nPoses + threads - 1) / threads;
        scoreKernel<<<blocks, threads>>>(dRec, dRecT, dRecR, nRec, dLig, dLigT, dLigR, nLig, dX,
                                         nPoses, dOut);
        ok = (cudaGetLastError() == cudaSuccess) && (cudaDeviceSynchronize() == cudaSuccess);
        if (ok)
            cudaMemcpy(outScores, dOut, nPoses * sizeof(float), cudaMemcpyDeviceToHost);
    }

    cudaFree(dRec);  cudaFree(dRecT); cudaFree(dRecR);
    cudaFree(dLig);  cudaFree(dLigT); cudaFree(dLigR);
    cudaFree(dX);    cudaFree(dOut);
    return ok;
}

}  // namespace stimlab::docking
