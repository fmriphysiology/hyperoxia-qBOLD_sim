% Script to prepare figures for sqBOLD/hqBOLD simulations

% Directory containing Monte Carlo results
% Assumed to be in same directory as 'figures/'
simdir='../storedprotonphases/'

%FIGURE 2

figure_qbold_nonoise(simdir)

%FIGURE 3

figure_qbold_noise_effect(simdir)

%FIGURE S3

figure_dbv_noise

