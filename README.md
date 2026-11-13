# Integrated Production Modelling

A Python workflow that links the reservoir, the wellbore and the simulator, starting with a nodal analysis engine (PVT, IPR, VLP) that reads its fluid straight from a simulation deck and is benchmarked / validated against PROSPER.

## Background

Nodal analysis shows where the reservoir (IPR) meets the tubing (VLP). Most of us run it in commercial software like PROSPER. This series asks: can the same logic be rebuilt in Python, with the fluid taken from the same deck the simulator uses, and still match PROSPER?

## Series Progress

* ✅ Part 1: Nodal Analysis Engine (PVT, IPR, VLP) benchmarked against PROSPER
* ⬜ Part 2: Python-generated VFP tables fed into ECLIPSE / OPM Flow
* ⬜ Later parts (planned, subject to change): scenario screening, gas lift optimisation, network modelling

## Part 1: What It Does

* **PVT** read directly from the ECLIPSE deck exported from Petrel (`PVDO`, `RSCONSTT`, `DENSITY`): Pb 80 bar, 26 API, reservoir temperature 76.85 C. Below Pb, Standing and Beggs-Robinson shapes are scaled to pass exactly through the deck values at Pb, with a bubble point and Bo that respond to temperature
* **IPR**: PI entry above Pb, Vogel below
* **VLP**: Beggs & Brill, marched from wellhead to bottomhole (vertical well, 8,858 ft, 3.958 in tubing)
* **Operating point** from the IPR/VLP intersection
* **Benchmark scripts** that read PROSPER's own exports (`.PRN`), so no numbers are typed by hand
* **3D GIFs**: VLP surface vs IPR sheet over liquid rate and water cut, and oil viscosity in pressure-temperature space

## Tools Used

* Python 3.10 (numpy, scipy, matplotlib, imageio; PyVista for visualization)
* PROSPER (benchmark)
* Petrel (deck export)

## Requirements

```
pip install numpy scipy matplotlib imageio pyvista
```

## About

Built as part of a self-directed learning project bridging reservoir simulation and production modelling. Documented step-by-step on LinkedIn; follow along for upcoming parts.

This is a work in progress, shared for learning purposes. Feedback and suggestions are very welcome!
