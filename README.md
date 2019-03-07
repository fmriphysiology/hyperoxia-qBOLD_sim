# Code developed to investigate systematic error in ASE based qBOLD measurements

This repository consists of two streams of work: (i) code to perform Monte Carlo simulations of the signal decay that occurs in the qBOLD experiment and (ii) code to generate figures based on the results of the former simulations. The results of the Monte Carlo simulations can be found here at https://doi.org/10.5287/bodleian:mvPY99a9D and is required to generate the figures. Code was developed using MATLAB version 8.4.0.150421 (R2014b). Please reference this code using the forthcoming Zenodo DOI, or for the dataset the DOI above, if you use either in your work.

## Monte Carlo simulations

The Monte Carlo approaches described by Boxerman et al. (1995) and Dickson et al. (2010) were implemented in MATLAB (simplevesselsim.m). This function is called by a simple script describing the physiology to be simulated e.g. simplevesselsim_jobR.m. The simulation code is also passed a structure containing physiological and physical properties of the system. 

## Figures

Each of the figure parts can be generated by the script prepare_figures.m. This script, and the others in the directory figures/ require the stored proton phases generated by the Monte Carlo simulation code (found here http://doi.org/XXXXX), which by default must be placed in the same directory as the figures/ directory. These phases can be appropriately combined to generate the signal decay of gradient echo (GRE), spin echo (SE), asymmetric spin echo (ASE) and gradient echo sampling of spin echo (GESSE) pulse sequences. 

## References

Boxerman, J.L., Bandettini, P.A., Kwong, K.K., Baker, J.R., Davis, T.L., Rosen, B.R., Weisskoff, R.M., 1995. The intravascular contribution to fMRI signal change: Monte Carlo modeling and diffusion-weighted studies in vivo. Magn. Reson. Med. 34, 4–10.

Dickson, J.D., Ash, T.W.J., Williams, G.B., Harding, S.G., Carpenter, T.A., Menon, D.K., Ansorge, R.E., 2010. Quantitative BOLD: the effect of diffusion. J. Magn. Reson. Imaging 32, 953–961. 
