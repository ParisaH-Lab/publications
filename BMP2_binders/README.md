
# From Code to Cure: BMP‑2 Binder Design Pipeline

This repository contains the computational workflows used in the study:

**From Code to Cure — Computationally Designed BMP‑2 Binders Using AI‑Integrated Pipelines for Controlled Bone Regeneration**  
Authors: Arash Asgari, Ben Burress‑Irving, et al.  
Corresponding: parisah@uoregon.edu

---

## Overview

We implemented a **two‑phase pipeline** to design de novo protein binders targeting the BMP‑2 knuckle epitope:

- **Phase I:** Rosetta/PyRosetta β‑strand motif grafting and docking.
- **Phase II:** Deep‑learning refinement (partial RFdiffusion + ProteinMPNN) and AlphaFold2 validation.

This repo includes **Phase I scripts** for motif grafting, sequence optimization, relaxation, and docking.

## Script Overview

- `dock_strands.py` — Python script for docking designed β-strand motifs onto the BMP2 scaffold, ensuring proper orientation and interface contacts.
- `minimize_strands.py` — Performs energy minimization on docked strand complexes to relieve clashes and optimize backbone geometry.
- `motif_graft.xml` — RosettaScripts protocol for grafting a functional motif onto a BMP2 binder scaffold using structural alignment and compatibility filters.
- `optimize_sequence2.xml` — RosettaScripts protocol for sequence optimization of the grafted binder, applying design movers and energy-based filters to improve binding affinity.


