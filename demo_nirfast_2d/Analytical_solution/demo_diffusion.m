addpath(genpath('C:\Users\letiz\OneDrive\Documents\Analytical solution'))
addpath(genpath('C:\Users\letiz\OneDrive\Documents\01.03.2022'))

% demo code: analytical solution to diffusion equations
% jingjing jiang 
% created: 2021.12.17 jingjing jiang
tstart = 0;
tstep = 1e-11;% unit: second
tend = 3e-9;
%time =tstart+tstep/2:tstep:tend;
time=0:tstep:tend;

% optical properties of the diffusive medium
mua = 0.01; % absorption [mm-1]
mus_p = 1; % reduced scattering [mm-1]
%c0=299792458000; % speed of light [mm/s] in vacuum
c0=3e11;        %->in set_mesh.m is like that, bsp. line 66
ri = 1.43;
v=c0/ri;        % speed of light in medium

srcpos = [-11 20 0]; % source positions [x y z] [mm]     % 2D or 3D, but chnage accordingly in tddiffusion, line 32
detpos = [-1.9502 20 0]; % detector positions [x y z] [mm]
Phi = tddiffusion(mua, ...
    mus_p, v, ...
    0.493, ...     
    srcpos, detpos+1, time);
%Letizia: put 0.493 as in the paper for the analytical solution (instead of 0)
                    %-> better matching with reconstructed signal

figure(1);
plot(time*1e9,Phi,'r');
xlabel('time [ns]') 
title('time of flight signal')