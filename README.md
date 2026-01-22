# Dirac-Equation Synchronization Dynamics (DESD)

This repo contains the MATLAB codes relevant to the paper titled "Designing topological cluster synchronization patterns with the Dirac operator" [1]:

* Dynamics.m: the main code required to run the DESD on a simple network. The inputs, outputs, and parameters are detailed in the comments at the top of the code.

* DESD.m: the auxiliary function called by Dynamics.m that implements the vector field associated with DESD.

* Spectra.m: plots the spectral cumulative density of the Hamiltonian operator defined for the user's choice of the network and mass constant. This can be used to visually determine the existence and number of isolated eigenstates that can induce topological cluster synchronization.

* DESD on SBM:
  - Plot_SBM.m: plots the data in DESD_m1_E1_SBM.mat that shows the topological cluster-synchronization induced by the DESD on an instance of the stochastic block model (SBM).
  - SBM: contains the code NetGen.m used to generate the SBM and the auxiliary function block_model.m

* DESD on Structural Connectome:
  - Plot_SC.m: plots the data in DESD_m1_E1_SC.mat and DESD_m1_E2_SC.mat that shows the topological cluster-synchronization induced by the DESD on the human structural connectome.
  - Data: contains the raw structural connectome data in SCmatrices88healthy.mat [2] and the preprocessing code Data_extraction.m

The codes can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed by the authors in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If you use these codes please cite

[1] A.A.A. Zaid and G. Bianconi, "Designing topological cluster synchronization patterns with the Dirac operator". Physical Review E (2026).
DOI: https://doi.org/10.1103/v65b-3jx7

If you use the data in the SCmatrices88healthy.mat file please cite:

[2] Antonín Škoch, Barbora Rehák Bučková, Jan Mareš, Jaroslav Tintěra, Pavel Sanda, Lucia Jajcay, Jiří Horáček, Filip Španiel, and Jaroslav Hlinka. Human
brain structural connectivity matrices–ready for modelling. Scientific Data, 9(1):486, 2022.

(c) A.A.A. Zaid (a.a.a.zaid@qmul.ac.uk) Ginestra Bianconi (g.bianconi@qmul.ac.uk)
