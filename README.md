# Ashkin_FDTD_and_Fluid_Solver

This code solves for deformation of a fluid-fluid interface using direct numerical simulation and front tracking. The CFD is implemented using the same method as 
in tryggvasson, who authored a good tutorial at: https://www3.nd.edu/~gtryggva/MultiphaseDNS/DNS-Solver.pdf and

https://pdfs.semanticscholar.org/f86a/c68eba685c5562837de229b3e6c1eaa8a9ef.pdf .

The electrodynamics is solved using FDTD in combination with direct differentiation of the electromagnetic stresss tensor. This was based on the publication I did at http://iopscience.iop.org/article/10.1088/0034-4885/78/12/122401/meta and a thesis at https://open.library.ubc.ca/media/download/pdf/24/1.0313408/3

another publication that used a similar code to this one is at  https://www.nature.com/articles/s41467-018-05706-3

all relevant documents should also be included in this folder.
