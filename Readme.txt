Read.me FLIPS

'Main' file contains the general example/problem of this package. The problem presented is a denoising problem with the Sparse Coding problem 
Min_f       ||f||_1 
s.t.        ||x – Df||_2 <= epsilon

The problem is solved with the oracles and state-of-the-art algorithms presented in the paper.
The other files contain the following information:
	- C_SALSA.m --> Contains the C-SALSA solver as described in:
	Manya V. Afonso, José M. Bioucas-Dias, and Mário A. T. Figueiredo. An Augmented
	Lagrangian Approach to the Constrained Optimization Formulation of Imaging Inverse
	Problems. IEEE Transactions On Image Processing, 2009
	- ChambollePock.m --> Contains the Chambolle-Pock solver as described in:
	Antonin Chambolle and Thomas Pock. On the ergodic convergence rates of a first-order
	primal–dual algorithm. Mathematical Programming, 159(1-2):253–287, 9 2016
	- DCT.m --> The file in which the square DCT-dictionary is computed. 
	- etafunc.m --> Contains the computations of the cost function eta, as described in the 								
	paper.
	- Fista Package --> Package from Tiep Vu https://github.com/tiepvupsu/FISTA open-source 	
	thv102@psu.edu, 4/6/2016.
	- FLIPS_Solver --> Solver FLIPS for the Quadratic Oracle
	- Frank_Wolf --> Solver FLIPS for the Linear Oracle
	- g_descendireciton_FW --> Part of the FLIPS Solver for Linear Oracle.
	- gradient_eta --> Contains the computations for the gradient of eta, as described in the 
	paper.
	- h_updatestep.m --> Update step of the variable h.
	- inputs --> Contains some standart images that can be used as input. The references of the inputs are listed in the report.
	- Main.m --> As described above.
	- patch2image.m --> Function that recreates the image from sliding image patches.
	- PGD_Oracle.m --> Solver for only Projected Gradient Descent.
 	- ProjectOntoL1Ball.m --> Projection function to ||.||_1 norm from:
	John Duchi, Shai@tti-C Org, and Tushar Chandra. Efficient Projections onto the l1-Ball
	for Learning in High Dimensions Google, Mountain View, CA 94043 Shai Shalev-Shwartz.
	Technical report, Proceedings of the 25th International Conference on Machine Learning,
	2008 -->https://stanford.edu/~jduchi/projects/DuchiShSiCh08.html
	- soft.m --> Soft thresholding function
	- stepsize_selection.m --> Exact line search function.

These files are allowed to be adjusted. However, without permission of the authors, it is not allowed to publish or distribute these files. 



