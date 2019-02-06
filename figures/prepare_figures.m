% Script to prepare figures for ASE qBOLD simulations

% Directory containing Monte Carlo results
% Assumed to be in same directory as 'figures/'
simdir='../storedprotonphases/'

%FIGURES S1, S2 & S3
figure_speedup_methods(simdir);

%FIGURES 2 & S4
figure_compare_gesse_ase(simdir);

%FIGURES 3, 4, 5, 6 & S5
figure_qbold_effects(simdir);

%FIGURES 7 & 8
figure_qbold_effects_sharan(simdir);

