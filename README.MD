## Py4siesta Utilities

### Based on NanoCore libraries (recent version: [**OORINano**](https://github.com/yhkimlab/OORINano) from YHKLab) 


This repository provides helper scripts for SIESTA calculations, wrapping and extending the NanoCore Python tools.

### 1. Set up the environment

Prepare the NanoCore library by moving the Python library directory to your working environment.

### 2. Prepare reference calculations

Create the directory `origin` containing:

- **`input/`** – SIESTA input (`.fdf`) files  
- **`slm_siesta_run`** – Slurm job script for running SIESTA

### 3. Run Py4siesta

Execute the main script:

```bash
python py4siesta.py

