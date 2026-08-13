// agent/SystemPrompt.h - the BioCAD assistant's system prompt.
//
// This is where the SAFETY BOUNDARY is encoded for the model: BioCAD predicts
// what a compound IS and DOES; it must never produce synthesis / route /
// precursor / manufacturability content. The UI layer appends the live panel
// list (so the model knows valid panel ids for the navigate/highlight tools).
#pragma once

#include <string>

namespace biocad::agent {

// The frozen safety + role preamble. AppShell appends the current panel catalog.
std::string safetySystemPrompt();

}  // namespace biocad::agent
