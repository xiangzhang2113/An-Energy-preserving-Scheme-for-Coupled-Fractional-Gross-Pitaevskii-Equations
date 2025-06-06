MATLAB code for the paper “An Energy-preserving Scheme for Coupled Fractional Gross-Pitaevskii Equations Based on Energy Discretization”.

###  Code Overview

- `Code_for_Tables_and_Figures.m`  
  Main script to reproduce all tables and figures in the paper.

- `fgp_solver.m`  
  Solves the Coupled Fractional Gross-Pitaevskii Equations and outputs the numerical solutions `u` and `v`.  
  This is the core solver used in all simulation scripts.

- `GP_ErrM_ErrE.m`  
  Computes the **mass and energy errors** for the Coupled Fractional Gross-Pitaevskii Equations under different numerical settings.

- `GP_M_E.m`  
  Computes the **mass**, **energy**, as well as **mass and energy errors** for the Coupled Fractional Gross-Pitaevskii Equations.  
  This script is used to validate the conservation properties of the proposed numerical scheme.

- `GP_Richardson.m`  
  Applies a **high-order explicit Richardson extrapolation method** to improve accuracy, and computes the resulting **mass**, **energy**, and their **errors**.  
  This file is used to demonstrate the conservation properties after applying explicit Richardson extrapolation.

