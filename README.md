# CoNbS 1Q/3Q Magnetic Structure Simulation Pipeline

This repository contains a Julia-based pipeline to simulate and analyze the magnetic structures (1Q and 3Q domains) of CoNbS using `Sunny.jl`. 

## Pipeline Overview

The computational pipeline consists of three main steps:
1. **1Q Domain Generation**: Deterministically generates the 1Q magnetic structure using the predefined sub-lattice relation `Spin(x+1, y+1, z, 2) = Spin(x, y, z, 1)`.
2. **3Q Domain Generation**: Minimizes a random spin configuration to stabilize a 3Q state, driven by biquadratic interactions.
3. **Visualization**: Takes the resulting `SampledCorrelations` data and visualizes the integrated structure factors along the Brillouin zone using `CairoMakie`.

## Scripts
- `src/CoreSimulation.jl`: A core Object-Oriented (OOP) module defining `SimulationParams`, `SystemConfig`, and `SimulationRunner`. It implements parallel thermal trajectories using `Threads.@threads` and reliably extracts `SampledCorrelations`.
- `calc_1Q_correlations.jl`: Instantiates the 1Q domain generation process and saves the output to a `.jld2` dataset.
- `calc_3Q_correlations.jl`: Instantiates the 3Q domain generation process and saves the output to a `.jld2` dataset.
- `plot_correlations.jl`: Processes the output datasets, applies mirror symmetry/Gaussian broadening, and exports figures.
- `run_pipeline.sh`: A SLURM bash script designed to allocate CPU cores and automatically run the pipeline in sequence.

## Environment & Run Details
- Multithreading: The workload is massively accelerated by parallelizing the collection of statistically independent samples.
- Usage: `sbatch run_pipeline.sh`
