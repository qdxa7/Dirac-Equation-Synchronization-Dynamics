# Dirac-Equation-Synchronization-Dynamics (DESD)

This repo contains the codes relevant to the paper titled "Designing topological cluster synchronization patterns with the Dirac operator".

DS_Dynamics1s.m: runs DESD and adiabatically increasing the value of the coupling strength sigma. Outputs the dynamical variables as a function of sigma. Uses random initial conditions.

DS_Dynamics1t.m: runs DESD for a fixed coupling strength sigma. Outputs the dynamical variables as a function of time. Uses random initial conditions.

DS_Dynamics2s.m: same as DS_Dynamics1s.m but uses an initial condition that is close to the selected eigenstate. This code is used for the stability analysis.

DS_Dynamics2t.m: same as DS_Dynamics1t.m but uses an initial condition that is close to the selected eigenstate. This code is used for the stability analysis.

The "Networks" subfolder contains the codes used to generate the Poisson (random) networks and the stochastic block models (SBM).

The "figure#" subfolders contain the codes used to plot the respective figures in the paper. They may also contain a shorter version of the DS_Dynamics*.m codes with only the lines required for the generation of the data used in the respective figures.

The code that generate the networks outputs a .mat file (starting with "RG" for a random network and "SBM" for an SBM) that contains the network data. The DS_Dynamics*.m codes take this .mat file as input and they output another .mat file (starting with "DS_Dyn") that contains the dynamics data. The figure plotting codes may take either .mat files as inputs depending on the plot. This is clarified by looking at the file name passed to the "load()" command in each code.

The codes can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed by the authors in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If you use these codes please cite

[1] A.A.A. Zaid and G. Bianconi, "Designing topological cluster synchronization patterns with the Dirac operator" (2025).

(c) A.A.A. Zaid (a.a.a.zaid@qmul.ac.uk) Ginestra Bianconi (g.bianconi@qmul.ac.uk)
