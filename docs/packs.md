# Data packs

The compound and target catalog is data, not code. A pack is one versioned JSON document
describing compounds and receptor targets; the application has no built-in knowledge of
CNS pharmacology or of any other domain - the domain is whatever the loaded packs say it
is (`src/packs/Pack.h:1-11`).

Packs are loaded by `biocad::packs` (`src/packs/Pack.h`, `src/packs/Pack.cpp`), cached
behind the docking preset accessors (`src/modules/docking/Presets.cpp:16-54`), and browsed
in the Presets panel (`src/ui/Panels.cpp:1477-1564`).

A pack describes what a compound **is** and where a ligand **binds**. There is
deliberately no synthesis, route, precursor or manufacturability field, and a loader that
encountered one would ignore it (`src/packs/Pack.h:9-11`).

## Schema version 1, annotated

`kSchemaVersion` is `1` (`src/packs/Pack.h:29`). Every field below is exactly what
`parse()` reads; anything else in the document is ignored.

```jsonc
{
  // REQUIRED, integer. Absent or non-integer -> parse error. Not equal to 1 -> unsupported
  // error. Never a silent skip.
  "schemaVersion": 1,

  // REQUIRED, non-empty string. The override key: a later pack with the same id replaces
  // this one wholesale.
  "id": "cns-monoamine",

  // REQUIRED, non-empty string. Shown in the Presets panel.
  "title": "CNS monoamine pharmacology",

  // Optional string.
  "description": "Reference compounds and receptor/enzyme targets ...",

  // Optional array. Absent is fine (safety-offtarget.json ships zero compounds).
  // Present but not an array -> parse error.
  "compounds": [
    {
      "id": "caffeine",                 // REQUIRED, non-empty. Stable slug; the join key.
      "name": "Caffeine",               // REQUIRED, non-empty. Display name.
      "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  // REQUIRED, non-empty.
      "drugClass": "Xanthine stimulant",         // Optional.
      "legalUs": "Unscheduled (US)",             // Optional; "legalStatus" is accepted as
                                                 // a fallback spelling.
      "notes": "Adenosine receptor antagonist; common reference stimulant.",  // Optional.

      // Optional object; every member optional. Non-object -> ignored entirely.
      // Later phases join on these; nothing in the loader dereferences them.
      "xrefs": {
        "chembl": "CHEMBL113",
        "pubchemCid": 2519,             // integer only; a quoted number is ignored
        "inchikey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
        "uniprot": ""                   // parsed for compounds too, but meaningful on targets
      },

      // Optional object of AUTHORED property values. Omit it and the consumer computes
      // the properties from the SMILES with the in-house chem engine instead of
      // inventing them. The block counts as authored only when molWeight > 0.
      "properties": {
        "formula": "C8H10N4O2",
        "molWeight": 194.19,            // > 0 is what sets hasProperties
        "logP": -0.07,
        "tpsa": 58.4,
        "hbd": 0,
        "hba": 6,
        "rotatableBonds": 0
      }
    }
  ],

  // Optional array. Present but not an array -> parse error.
  "targets": [
    {
      "id": "DAT",                                  // REQUIRED, non-empty.
      "name": "DAT (dopamine transporter)",         // REQUIRED, non-empty.
      "pdb": "4M48",                                // Optional string - but REQUIRED if
                                                    // "box" is present. null/absent means
                                                    // "no structure".
      "uniprot": "Q01959",                          // Optional.
      "chembl": "CHEMBL238",                        // Optional.
      "headline": true,                             // Optional, default false. Headline
                                                    // receptors are prepared during
                                                    // provisioning.

      // Optional object. Present -> this target is dockable. Centre defaults to the
      // origin; each edge defaults to 22.0 A; any edge <= 0 is a parse error.
      "box": { "cx": -2.1, "cy": 12.4, "cz": -6.8, "sx": 20.0, "sy": 20.0, "sz": 20.0 },

      // Optional array of strings; non-strings inside are skipped. Panel membership tags,
      // e.g. an off-target screening panel.
      "panels": ["safetyscreen44"]
    }
  ]
}
```

## Required fields and what happens when they are missing

Every required field goes through `requireString` (`src/packs/Pack.cpp:33-39`), which
treats missing, non-string and empty-string identically and throws
`Error::parse("<sourcePath>: missing required string field \"<key>\"")`. `Error` is the
named exception type from `src/core/Error.h`.

|Field|Required by|Failure|
|---|---|---|
|`schemaVersion`|`parse`, `src/packs/Pack.cpp:148-151`|`Error::parse`, `missing integer "schemaVersion"`|
|`schemaVersion == 1`|`src/packs/Pack.cpp:154-158`|`Error::unsupported`, naming the version found and the version this build reads|
|top-level `id`|`src/packs/Pack.cpp:160`|`Error::parse`|
|top-level `title`|`src/packs/Pack.cpp:161`|`Error::parse`|
|compound `id`|`src/packs/Pack.cpp:62`|`Error::parse`|
|compound `name`|`src/packs/Pack.cpp:63`|`Error::parse`|
|compound `smiles`|`src/packs/Pack.cpp:64`|`Error::parse`, with the compound id in the message - a compound without SMILES cannot be computed on|
|target `id`|`src/packs/Pack.cpp:87`|`Error::parse`|
|target `name`|`src/packs/Pack.cpp:88`|`Error::parse`|
|`pdb`, when `box` is present|`src/packs/Pack.cpp:110-113`|`Error::parse`, "a box without a PDB id is an unverifiable binding site"|

Structural failures throw the same way: a non-object document
(`src/packs/Pack.cpp:146`), a non-array `compounds` (`:167`) or `targets` (`:181`), and
text that is not JSON at all (`parseString`, `:196-198`). `tests/test_packs.cpp:112-119`
covers the required-field set.

## The hard rules

Each rule below is implemented at the cited line, not merely intended.

|Rule|Implementation|
|---|---|
|An unknown `schemaVersion` is a load error, surfaced in the Presets panel, never a silent skip|`src/packs/Pack.cpp:154-158` throws `Error::unsupported`; `loadFrom` catches it into `LoadReport::errors` (`:226-230`); the panel renders `errors` in red before anything else (`src/ui/Panels.cpp:1492-1500`). Rationale at `src/packs/Pack.h:26-28`: a pack that silently vanished is indistinguishable from a broken application.|
|A pack with no `schemaVersion` at all is equally a hard error|`src/packs/Pack.cpp:148-151`; `tests/test_packs.cpp:84-86`|
|A `box` requires a `pdb`|`src/packs/Pack.cpp:110-113`; `tests/test_packs.cpp:89-95`|
|A degenerate box is rejected, not clamped into something plausible|`src/packs/Pack.cpp:104-106` - any of `sx`, `sy`, `sz` `<= 0` throws; `tests/test_packs.cpp:104-109`|
|A target without a box is an honest coverage gap and is not dockable|`hasBox` stays false (`src/packs/Pack.cpp:96-109`, `src/packs/Pack.h:68-73`); the preset cache skips it with `if (!t.hasBox) continue;` (`src/modules/docking/Presets.cpp:28`), so it never reaches `targetPresets()`, `findPreset()` or `presetNames()`; the panel draws it grey as "no box" (`src/ui/Panels.cpp:1549-1551`). `tests/test_packs.cpp:97-102`|
|Duplicate compound ids inside one pack are rejected at parse|`src/packs/Pack.cpp:171-174`; `tests/test_packs.cpp:121-127`|
|Duplicate target ids inside one pack are rejected at parse|`src/packs/Pack.cpp:185-188`|
|Duplicate ids **across** packs are reported in `LoadReport::errors` and the later one is dropped|Compounds `src/packs/Pack.cpp:295-304`, targets `:314-323`. The message names the id and both pack ids. This is why `compounds()` and `targets()` are non-const (`src/packs/Pack.h:98-103`). `tests/test_packs.cpp:148-159`|
|A whole pack is overridden by pack id, later-wins|`effectivePacks` (`src/packs/Pack.cpp:273-286`) keeps first-seen position and replaces the pointer. The replacement is wholesale, not a merge, so the user can always name the one document that produced a row (`:270-272`). Overriding is **not** an error. `tests/test_packs.cpp:129-146`|
|A missing pack directory is not an error|`loadFrom` returns an empty report when the directory does not exist (`src/packs/Pack.cpp:203-204`); `tests/test_packs.cpp:171-175`|
|Having **no** packs at all is an error|`loadBuiltin` synthesises the message "no built-in packs found - the application has no compound or target catalog" (`src/packs/Pack.cpp:256-260`)|
|An unreadable file is reported per-file and does not abort the load|`src/packs/Pack.cpp:216-219`|

`findTarget` resolves case-insensitively by id, then by display name, then by name prefix,
so `"DAT"` and `"DAT (dopamine transporter)"` both resolve
(`src/packs/Pack.cpp:328-343`; `tests/test_packs.cpp:161-169`). The same three-step
resolution exists over the dockable subset in `findPreset`
(`src/modules/docking/Presets.cpp:56-67`).

## Resolution order

`loadBuiltin()` (`src/packs/Pack.cpp:254-266`):

1. **Built-in packs**, `assets/packs/*.json` beside the executable. `builtinPackDir()`
   (`:235-252`) looks for `<exeDir>/assets/packs`, then the dev-tree fallback three
   directories up from the exe (`build/<preset>/bin/BioCAD.exe`), then a relative
   `assets/packs`. These are flagged `builtin = true`.
2. **User packs**, `AppPaths::instance().root() / "packs"` - that is
   `%APPDATA%\BioCAD\packs\*.json`. Flagged `builtin = false`.

Within each directory, `loadFrom` collects every `*.json` regular file and **sorts by
filename** before parsing, so the load order is deterministic and independent of the
filesystem's enumeration order (`src/packs/Pack.cpp:206-212`). Errors from both
directories are concatenated into one `LoadReport::errors`.

Because user packs are appended after built-ins and `effectivePacks` is later-wins, a user
pack whose `id` matches a built-in replaces it entirely. A user pack with a fresh `id` is
additive.

## Author your own pack

1. Create `%APPDATA%\BioCAD\packs` if it does not exist.
2. Write a `.json` file there - the filename only affects load order, the `id` field is
   what identifies the pack. Pick an `id` that does not collide with a shipped pack unless
   you intend to replace it wholesale.
3. Start from the minimum viable document:

   ```json
   {
     "schemaVersion": 1,
     "id": "my-lab",
     "title": "My lab compounds",
     "compounds": [
       { "id": "my-cmpd-1", "name": "Compound 1", "smiles": "CC(N)Cc1ccccc1" }
     ]
   }
   ```

   With no `properties` block, molecular weight, logP, TPSA and the rest are computed from
   the SMILES by the in-house chem engine rather than being authored
   (`src/packs/Pack.h:52-53`, `src/packs/Pack.cpp:133-141`).
4. Add a target only if you can justify its binding site. With a `pdb` and a `box` it
   becomes dockable; with a `pdb` and no `box`, or with neither, it is listed as a
   coverage gap:

   ```json
   "targets": [
     { "id": "MYT", "name": "My target", "pdb": "1ABC",
       "box": { "cx": 10.2, "cy": -4.0, "cz": 3.1, "sx": 22.0, "sy": 22.0, "sz": 22.0 } }
   ]
   ```

   Omitting `pdb` while supplying `box` is rejected. Box edges default to 22.0 A each if
   you omit them, and must be positive.
5. Set `"headline": true` on at most a few targets: headline receptors are prepared during
   provisioning rather than on first use (`src/packs/Pack.h:74`,
   `src/modules/docking/Presets.cpp:30`).
6. Open the **Presets** panel and press **Reload packs**
   (`src/ui/Panels.cpp:1487`, `docking::reloadTargetPacks()` at
   `src/modules/docking/Presets.cpp:50-54`). This re-reads every pack from disk without
   restarting the application.
7. Read the panel top-down: the red **LOAD ERRORS** block first, then **LOADED PACKS**
   (id, built-in vs user, compound count, target count, source path), then **TARGETS**
   with each target's PDB, `dockable`/`no box` badge, and tags. If your pack is not in the
   second table, its reason is in the first.

## Shipped packs

Counts read from `assets/packs/*.json` on disk.

|File|`id`|Title|Compounds|Targets|Dockable (boxed)|Headline|
|---|---|---|---:|---:|---:|---:|
|`analgesics-otc.json`|`analgesics-otc`|Analgesics, antipyretics and common OTC drugs|18|4|0|0|
|`cns-monoamine.json`|`cns-monoamine`|CNS monoamine pharmacology|31|29|29|4|
|`metabolic-and-nutrition.json`|`metabolic-and-nutrition`|Metabolic, nutritional and performance compounds|18|4|0|0|
|`safety-offtarget.json`|`safety-offtarget`|Safety and off-target panel|0|22|0|0|
|**Total**|||**67**|**59**|**29**|**4**|

The shape of that table is the point: 59 targets are catalogued and 29 of them are
dockable. The other 30 have no verified binding-site box, so they are listed as coverage
gaps rather than handed a fabricated site. `safety-offtarget.json` is the extreme case -
22 targets, none boxed, zero compounds - which is an honest statement that the off-target
panel is enumerated but not yet dockable, and is exactly the kind of fact the Presets
panel exists to make visible.
