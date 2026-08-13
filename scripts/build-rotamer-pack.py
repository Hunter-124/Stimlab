#!/usr/bin/env python3
"""Derive the backbone-dependent rotamer pack in assets/packs/rotamers/ from real
crystal structures.

Why this script exists at all: the Dunbrack 2010 smoothed backbone-dependent
rotamer library is CC BY 4.0, but the only distribution channel is a download page
released after a manually approved licence application
(https://dunbrack.fccc.edu/lab/bbdep2010_download_page_example). There is no
unauthenticated URL, so the library cannot be vendored here. Inventing rotamer
angles instead is not an option, so the pack is DERIVED, with the derivation
recorded in the pack itself: every mean angle, bond length, bond angle and
probability below is measured from the PDB entries listed in `dataset.entries`.

The binning convention (rotameric chi thirds around -60/180/+60, a non-rotameric
final chi for ASN/ASP/GLN/GLU/HIS/PHE/TYR/TRP, phi/psi bins) follows Dunbrack's;
the NUMBERS are ours and are not interchangeable with theirs.

Usage:  python3 scripts/build-rotamer-pack.py <dir-of-pdb.gz> <out.json>
"""

import gzip
import json
import math
import os
import sys
from collections import defaultdict

# --- chemistry tables (connectivity and dihedral definitions are chemical fact) ---

CHI = {
    "ARG": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"),
            ("CB", "CG", "CD", "NE"), ("CG", "CD", "NE", "CZ")],
    "ASN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "ASP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "CYS": [("N", "CA", "CB", "SG")],
    "GLN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "GLU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "HIS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")],
    "ILE": [("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")],
    "LEU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "LYS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"),
            ("CB", "CG", "CD", "CE"), ("CG", "CD", "CE", "NZ")],
    "MET": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD"), ("CB", "CG", "SD", "CE")],
    "PHE": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "SER": [("N", "CA", "CB", "OG")],
    "THR": [("N", "CA", "CB", "OG1")],
    "TRP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "TYR": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "VAL": [("N", "CA", "CB", "CG1")],
}

# Which trailing chi is non-rotameric (a continuous distribution, not three wells).
NON_ROTAMERIC = {"ASN": 2, "ASP": 2, "GLN": 3, "GLU": 3,
                 "HIS": 2, "PHE": 2, "TRP": 2, "TYR": 2}
# Symmetric terminal chi: the two branch atoms are indistinguishable, so chi and
# chi+180 are the same physical conformation and are folded into (-90, 90].
SYMMETRIC = {"ASP": 2, "PHE": 2, "TYR": 2, "GLU": 3}

# Side-chain build order: atom, its three already-placed parents (a, b, c) for the
# NeRF placement a-b-c->d, and the torsion source ("chi N" or "fixed").
BUILD = {
    "ALA": [("CB", ("N", "C", "CA"), 0)],
    "SER": [("CB", ("N", "C", "CA"), 0), ("OG", ("N", "CA", "CB"), 1)],
    "CYS": [("CB", ("N", "C", "CA"), 0), ("SG", ("N", "CA", "CB"), 1)],
    "THR": [("CB", ("N", "C", "CA"), 0), ("OG1", ("N", "CA", "CB"), 1),
            ("CG2", ("N", "CA", "CB"), 0)],
    "VAL": [("CB", ("N", "C", "CA"), 0), ("CG1", ("N", "CA", "CB"), 1),
            ("CG2", ("N", "CA", "CB"), 0)],
    "ILE": [("CB", ("N", "C", "CA"), 0), ("CG1", ("N", "CA", "CB"), 1),
            ("CG2", ("N", "CA", "CB"), 0), ("CD1", ("CA", "CB", "CG1"), 2)],
    "LEU": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("CD1", ("CA", "CB", "CG"), 2), ("CD2", ("CA", "CB", "CG"), 0)],
    "ASP": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("OD1", ("CA", "CB", "CG"), 2), ("OD2", ("CA", "CB", "CG"), 0)],
    "ASN": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("OD1", ("CA", "CB", "CG"), 2), ("ND2", ("CA", "CB", "CG"), 0)],
    "GLU": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("CD", ("CA", "CB", "CG"), 2), ("OE1", ("CB", "CG", "CD"), 3),
            ("OE2", ("CB", "CG", "CD"), 0)],
    "GLN": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("CD", ("CA", "CB", "CG"), 2), ("OE1", ("CB", "CG", "CD"), 3),
            ("NE2", ("CB", "CG", "CD"), 0)],
    "MET": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("SD", ("CA", "CB", "CG"), 2), ("CE", ("CB", "CG", "SD"), 3)],
    "LYS": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("CD", ("CA", "CB", "CG"), 2), ("CE", ("CB", "CG", "CD"), 3),
            ("NZ", ("CG", "CD", "CE"), 4)],
    "ARG": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("CD", ("CA", "CB", "CG"), 2), ("NE", ("CB", "CG", "CD"), 3),
            ("CZ", ("CG", "CD", "NE"), 4), ("NH1", ("CD", "NE", "CZ"), 0),
            ("NH2", ("CD", "NE", "CZ"), 0)],
    "HIS": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("ND1", ("CA", "CB", "CG"), 2), ("CD2", ("CA", "CB", "CG"), 0),
            ("CE1", ("CB", "CG", "ND1"), 0), ("NE2", ("CB", "CG", "CD2"), 0)],
    "PHE": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("CD1", ("CA", "CB", "CG"), 2), ("CD2", ("CA", "CB", "CG"), 0),
            ("CE1", ("CB", "CG", "CD1"), 0), ("CE2", ("CB", "CG", "CD2"), 0),
            ("CZ", ("CG", "CD1", "CE1"), 0)],
    "TYR": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("CD1", ("CA", "CB", "CG"), 2), ("CD2", ("CA", "CB", "CG"), 0),
            ("CE1", ("CB", "CG", "CD1"), 0), ("CE2", ("CB", "CG", "CD2"), 0),
            ("CZ", ("CG", "CD1", "CE1"), 0), ("OH", ("CD1", "CE1", "CZ"), 0)],
    "TRP": [("CB", ("N", "C", "CA"), 0), ("CG", ("N", "CA", "CB"), 1),
            ("CD1", ("CA", "CB", "CG"), 2), ("CD2", ("CA", "CB", "CG"), 0),
            ("NE1", ("CB", "CG", "CD1"), 0), ("CE2", ("CB", "CG", "CD2"), 0),
            ("CE3", ("CB", "CG", "CD2"), 0), ("CZ2", ("CG", "CD2", "CE2"), 0),
            ("CZ3", ("CG", "CD2", "CE3"), 0), ("CH2", ("CD2", "CE2", "CZ2"), 0)],
}

BIN = 60          # phi/psi bin width, degrees
MIN_BIN_COUNT = 20   # below this a phi/psi bin is not emitted; the caller falls back

# --- geometry ---


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(a):
    return math.sqrt(dot(a, a))


def dihedral(p0, p1, p2, p3):
    # IUPAC sign convention (the "praxeolitic" projection form): looking down
    # p1->p2, a clockwise rotation of p3 relative to p0 is positive. Getting this
    # sign backwards silently mirrors every rotamer (LEU would read pt instead of
    # mt), so it is fixed here once rather than per call site.
    b0 = sub(p0, p1)
    b1 = sub(p2, p1)
    b2 = sub(p3, p2)
    n = norm(b1)
    u = (b1[0] / n, b1[1] / n, b1[2] / n)
    v = tuple(b0[i] - dot(b0, u) * u[i] for i in range(3))
    w = tuple(b2[i] - dot(b2, u) * u[i] for i in range(3))
    return math.degrees(math.atan2(dot(cross(u, v), w), dot(v, w)))


def angle(p0, p1, p2):
    a, b = sub(p0, p1), sub(p2, p1)
    c = dot(a, b) / (norm(a) * norm(b))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def circ_mean(vals):
    s = sum(math.sin(math.radians(v)) for v in vals)
    c = sum(math.cos(math.radians(v)) for v in vals)
    return math.degrees(math.atan2(s / len(vals), c / len(vals)))


def circ_sd(vals):
    s = sum(math.sin(math.radians(v)) for v in vals) / len(vals)
    c = sum(math.cos(math.radians(v)) for v in vals) / len(vals)
    r = min(1.0, math.sqrt(s * s + c * c))
    if r < 1e-12:
        return 180.0
    return math.degrees(math.sqrt(-2.0 * math.log(r)))


def wrap(a):
    while a <= -180.0:
        a += 360.0
    while a > 180.0:
        a -= 360.0
    return a


def fold_symmetric(a):
    """Fold a two-fold symmetric torsion into (-90, 90]."""
    a = wrap(a)
    if a > 90.0:
        a -= 180.0
    elif a <= -90.0:
        a += 180.0
    return a


def bond_to_chi(resn):
    """Map the rotatable bond (parents[1], parents[2]) onto its chi index.

    Every atom placed about the same bond as chi_k rotates with chi_k, so its
    torsion is stored as an offset from chi_k rather than as a standalone mean.
    """
    out = {}
    for k, quad in enumerate(CHI.get(resn, []), start=1):
        out[(quad[1], quad[2])] = k
    return out


def well(a):
    """Rotameric third: p ~ +60, t ~ 180, m ~ -60."""
    a = wrap(a)
    if 0.0 < a <= 120.0:
        return "p"
    if a > 120.0 or a <= -120.0:
        return "t"
    return "m"


# --- PDB parsing (only what this script needs: one model, altLoc A/blank) ---


def read_pdb(path):
    chains = defaultdict(list)   # chain -> list of (resname, resid, {atom: xyz})
    cur = {}
    key = None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM  "):
                continue
            alt = line[16]
            if alt not in (" ", "A"):
                continue
            name = line[12:16].strip()
            resn = line[17:20].strip()
            ch = line[21]
            rid = line[22:27]
            occ = float(line[54:60]) if line[54:60].strip() else 1.0
            if occ < 0.5:
                continue
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            k = (ch, rid)
            if k != key:
                if key is not None:
                    chains[key[0]].append((curn, key[1], cur))
                key, cur, curn = k, {}, resn
            cur[name] = xyz
    if key is not None:
        chains[key[0]].append((curn, key[1], cur))
    return chains


def main():
    pdbdir, out = sys.argv[1], sys.argv[2]
    files = sorted(f for f in os.listdir(pdbdir) if f.endswith(".pdb.gz"))

    # residue -> bin key -> rotamer key -> list of chi tuples
    obs = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    # residue -> atom -> lists of (bond, angle, torsion-delta)
    geom = defaultdict(lambda: defaultdict(lambda: ([], [], [])))
    used_entries, nres = [], 0

    for f in files:
        chains = read_pdb(os.path.join(pdbdir, f))
        touched = False
        for ch, residues in chains.items():
            for i, (resn, rid, at) in enumerate(residues):
                if resn not in CHI and resn != "ALA":
                    continue
                if i == 0 or i + 1 >= len(residues):
                    continue
                prev, nxt = residues[i - 1], residues[i + 1]
                need = ("N", "CA", "C")
                if any(a not in at for a in need):
                    continue
                if "C" not in prev[2] or "N" not in nxt[2]:
                    continue
                phi = dihedral(prev[2]["C"], at["N"], at["CA"], at["C"])
                psi = dihedral(at["N"], at["CA"], at["C"], nxt[2]["N"])

                build = BUILD.get(resn)
                if not build:
                    continue
                if any(a not in at for a, _, _ in build):
                    continue

                chis, chis_raw = [], []
                if resn != "ALA":
                    ok = True
                    for j, quad in enumerate(CHI[resn], start=1):
                        if any(a not in at for a in quad):
                            ok = False
                            break
                        v = dihedral(*[at[a] for a in quad])
                        chis_raw.append(v)
                        if SYMMETRIC.get(resn) == j:
                            v = fold_symmetric(v)
                        chis.append(v)
                    if not ok:
                        continue

                # Internal coordinates for the build tree. An atom whose torsion is
                # about the SAME bond as a chi rotates with that chi, so what is
                # measured for it is the offset from that chi (LEU CD2 sits at
                # chi2 + ~120), not a standalone mean - a standalone mean over a
                # rotating branch atom is an average of every rotamer at once and
                # places the atom nowhere. The offset is taken against the RAW chi,
                # not the symmetry-folded one, or ASP OD2 would come out bimodal at
                # 0/180 and average onto OD1.
                bond_chi = bond_to_chi(resn)
                for atom, (pa, pb_, pc), _src in build:
                    if pa not in at or pb_ not in at or pc not in at:
                        continue
                    b, g, t = geom[resn][atom]
                    b.append(norm(sub(at[atom], at[pc])))
                    g.append(angle(at[pb_], at[pc], at[atom]))
                    tor = dihedral(at[pa], at[pb_], at[pc], at[atom])
                    k = bond_chi.get((pb_, pc), 0)
                    t.append(wrap(tor - chis_raw[k - 1]) if k else tor)

                if resn == "ALA":
                    continue
                nrot = NON_ROTAMERIC.get(resn, len(chis) + 1)
                rot = tuple(well(c) for j, c in enumerate(chis, start=1) if j != nrot)
                pb = int(math.floor((phi + 180.0) / BIN)) * BIN - 180
                sb = int(math.floor((psi + 180.0) / BIN)) * BIN - 180
                obs[resn][(pb, sb)][rot].append(tuple(chis))
                obs[resn][("any", "any")][rot].append(tuple(chis))
                nres += 1
                touched = True
        if touched:
            used_entries.append(f.split(".")[0].upper())

    residues = {}
    for resn in sorted(set(list(BUILD.keys()))):
        entry = {
            "chiCount": len(CHI.get(resn, [])),
            "nonRotamericChi": NON_ROTAMERIC.get(resn, 0),
            "symmetricChi": SYMMETRIC.get(resn, 0),
            "build": [],
            "backboneIndependent": [],
            "bins": [],
        }
        bond_chi = bond_to_chi(resn)
        for atom, parents, _src in BUILD[resn]:
            b, g, t = geom[resn][atom]
            if not b:
                continue
            k = bond_chi.get((parents[1], parents[2]), 0)
            entry["build"].append({
                "atom": atom,
                "parents": list(parents),
                "bondLength": round(sum(b) / len(b), 4),
                "bondAngle": round(sum(g) / len(g), 3),
                # When chi > 0 this is an OFFSET added to chi_k (0 for the atom that
                # defines chi_k itself); when chi == 0 it is the torsion outright.
                "torsion": round(circ_mean(t), 3),
                "torsionSd": round(circ_sd(t), 3),
                "chi": k,
                "observations": len(b),
            })

        def rotamers_for(bucket):
            total = sum(len(v) for v in bucket.values())
            outr = []
            for rot, rows in sorted(bucket.items(), key=lambda kv: -len(kv[1])):
                if len(rows) < 3:
                    continue
                nchi = len(rows[0])
                means, sds = [], []
                for j in range(nchi):
                    vals = [r[j] for r in rows]
                    means.append(round(circ_mean(vals), 2))
                    sds.append(round(circ_sd(vals), 2))
                outr.append({
                    "name": "".join(rot) if rot else "-",
                    "count": len(rows),
                    "probability": round(len(rows) / total, 6),
                    "chi": means,
                    "chiSd": sds,
                })
            return outr, total

        if resn in CHI:
            bi, bitotal = rotamers_for(obs[resn][("any", "any")])
            entry["backboneIndependent"] = bi
            entry["backboneIndependentCount"] = bitotal
            for (pb, sb), bucket in sorted(obs[resn].items(), key=lambda kv: str(kv[0])):
                if pb == "any":
                    continue
                rots, total = rotamers_for(bucket)
                if total < MIN_BIN_COUNT or not rots:
                    continue
                entry["bins"].append({"phi": pb, "psi": sb, "count": total, "rotamers": rots})
        residues[resn] = entry

    pack = {
        "schemaVersion": 1,
        "id": "rotamers-pdb-derived-2026",
        "libraryName": "BioCAD PDB-derived backbone-dependent rotamer library (2026-08)",
        "title": "Backbone-dependent rotamer means, probabilities and side-chain internal "
                 "coordinates measured from high-resolution crystal structures",
        "method": ("Every number in this file is measured. Side chains were taken from "
                   "X-ray entries at <= 1.4 A resolution, occupancy >= 0.5, altLoc blank or "
                   "A. chi angles use the IUPAC atom quadruples; rotameric chi are binned "
                   "into thirds (p ~ +60, t ~ 180, m ~ -60) and the terminal non-rotameric "
                   "chi of ASN/ASP/GLN/GLU/HIS/PHE/TRP/TYR is left continuous and reported "
                   "as a circular mean with its circular standard deviation. Symmetric "
                   "terminal chi (ASP chi2, GLU chi3, PHE/TYR chi2) are folded into "
                   "(-90, 90]. phi/psi are binned at " + str(BIN) + " degrees; a bin with "
                   "fewer than " + str(MIN_BIN_COUNT) + " observations is omitted and the "
                   "backbone-independent aggregate is used instead, which the consumer "
                   "reports as an assumption. Bond lengths and bond angles are the dataset "
                   "means for the same atoms."),
        "notDunbrack": ("This is NOT the Dunbrack 2010 smoothed backbone-dependent library. "
                        "That library is CC BY 4.0 but is distributed only through a "
                        "manually approved licence application, so it cannot be vendored "
                        "here. The binning convention below follows Shapovalov & Dunbrack "
                        "2011, Structure 19:844-858; the numbers do not come from it and "
                        "are not interchangeable with it."),
        "attribution": [
            "Coordinates: RCSB PDB (wwPDB), CC0 1.0. Entry ids in dataset.entries.",
            "Binning convention after Shapovalov & Dunbrack 2011, Structure 19:844-858.",
        ],
        "phiPsiBinDegrees": BIN,
        "minimumBinObservations": MIN_BIN_COUNT,
        "dataset": {
            "source": "RCSB PDB search: X-ray, resolution <= 1.4 A, 150-600 polymer residues",
            "entryCount": len(used_entries),
            "residuesMeasured": nres,
            "entries": used_entries,
        },
        "residues": residues,
    }
    with open(out, "w") as fh:
        json.dump(pack, fh, indent=1, sort_keys=False)
    print("entries", len(used_entries), "residues", nres,
          "bytes", os.path.getsize(out))


if __name__ == "__main__":
    main()
