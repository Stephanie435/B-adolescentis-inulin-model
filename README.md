# Inulin-*Bifidobacterium adolescentis* Consumer-Resource Model

## Project overview
This repository contains the MATLAB code used to simulate the effect of initial inulin availability on the population dynamics of *Bifidobacterium adolescentis*.

The model is a simplified consumer-resource system adapted from Jansma et al. (2026). It includes one bacterial population and one dietary fibre resource (inulin), with Monod-type saturation in the bacterial growth term.

## Model
The model is defined by:

$$
g(F)=\frac{c_0cF}{1+F/K_s}
$$

$$
\frac{dB}{dt}=B(g(F)-d)(k-B)
$$

$$
\frac{dF}{dt}=-cFB
$$

where:
- `B(t)` = scaled *B. adolescentis* abundance
- `F(t)` = relative inulin availability
- `c0` = resource-to-growth conversion coefficient
- `c` = interaction/consumption coefficient
- `Ks` = Monod half-saturation parameter
- `d` = bacterial loss parameter
- `k` = carrying capacity

The system is closed, so no additional inulin enters after the initial condition `F(0)`.

## Simulated inulin scenarios
Three initial inulin conditions are compared:
- Low: `F0 = 1.5`
- Baseline: `F0 = 3.0`
- High: `F0 = 6.0`

These are scaled relative resource units and do not correspond directly to grams of dietary inulin.

## MATLAB file
Run:

`Stephanie_Henning_27944611.m`

The script:
- solves the coupled ODE system using `ode45`;
- plots bacterial abundance over time;
- plots inulin depletion over time;
- compares peak bacterial abundance;
- performs a sensitivity analysis for `c = {0.70, 0.90, 1.00}` and `Ks = {1, 2, 4}`.

## Requirements
- MATLAB
- No additional toolboxes are required for the core simulations.

## Reproducing the results
1. Download or clone this repository.
2. Open `Stephanie_Henning_27944611.m` in MATLAB.
3. Run the script.
4. The figures and summary output will be generated automatically.




