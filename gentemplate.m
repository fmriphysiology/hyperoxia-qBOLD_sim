function p=gentemplate;

    p.R = 20e-6; %m
    p.D = 1.3e-9; %m^2/s
    p.vesselDensity = 0.05;
    p.B0 = 3; %T
    p.gamma = 2*pi*42.58e6;
    p.TE = 60e-3;
    p.deltaTE = 2e-3; 
    p.dt = 200e-6;
    p.deltaChi = 0.264e-6.*0.4.*(1-0.6);
    p.N = 10000;

return;