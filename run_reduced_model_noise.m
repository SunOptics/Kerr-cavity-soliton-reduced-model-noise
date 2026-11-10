%% run_reduced_model_noise.m
% Noise spectra of a stationary Kerr cavity soliton from the four-coordinate
% reduced model (amplitude eta, frequency shift Omega, position tau0,
% phase phi) of Y. Sun et al., arXiv:2605.14614 (2026).
%
% Edit the parameters in sections 1-4 and run the script. It calls
%   reduced_model_noise_solver.m  (all physics, no plotting)
%   plot_reduced_model_noise.m    (figure in the layout of Fig. 4 of the paper)
% and writes the figure (.png, .pdf) and the results (.mat) to ./output/.
%
% The default values are the fiber-ring example of Fig. 4(a),(b) of the
% paper (normalized P = 3, delta = 9), with quantum-limited noise.
%
% Author:    Yifan Sun (yifan.sun@ulb.be)
%            Service OPERA-Photonics, Université libre de Bruxelles (ULB),
%            B-1050 Brussels, Belgium
% Reference: Y. Sun, T. Bunel, S. Glazyrina, G. Semaan, F. Bretenaker,
%            S. Coen, S.-P. Gorza, and F. Leo, "Collective-Coordinate
%            Fluctuations of Driven-Dissipative Solitons",
%            arXiv:2605.14614 [physics.optics] (2026),
%            https://arxiv.org/abs/2605.14614
% Figures:   defaults -> Fig. 4(a).1-4 and 4(b).1-4 (P = 3, delta = 9);
%            P = 4, delta = 9  -> Fig. 4(c).1-4;
%            P = 4, delta = 19 -> Fig. 4(d).1-4
%            (set par.pump_power_unit = 'normalized' to enter P directly).

clear; close all;

%% 1) Resonator
par.L_m            = 380;        % roundtrip fiber length [m]
par.beta1_s_per_m  = 4.8823e-9;  % inverse group velocity n_g/c [s/m]; FSR = 1/(beta1*L)
par.beta2_s2_per_m = -20e-27;    % group-velocity dispersion [s^2/m] (-20 ps^2/km); must be < 0
par.gamma_per_W_m  = 1.8e-3;     % Kerr nonlinear coefficient [1/(W m)]
par.finesse        = 24;         % loaded finesse F = 2*pi/(roundtrip loss)
par.theta          = 0.1;        % power coupling ratio of the input coupler
par.tauR_s         = 3e-15;      % Raman time constant [s] (about 3 fs in silica)
par.lambda0_m      = 1550e-9;    % pump wavelength [m]

%% 2) Drive
par.pump_power      = 0.29512;   % pump power at the resonator input
par.pump_power_unit = 'W';       % 'W', or 'normalized' to give P of the paper directly
par.detuning        = 9;         % pump-resonance detuning, > 0 on the soliton side
par.detuning_unit   = 'normalized';  % 'normalized' (delta of the paper),
                                     % 'rad' (roundtrip phase) or 'Hz' (frequency)

%% 3) Noise
% Diffusion matrix D = noise_enhancement_ratio * D_vacuum. 1 = quantum
% limited; e.g. 10 injects ten times the quantum-limited noise into every
% coordinate. The analytic PSDs scale linearly with this ratio.
par.noise_enhancement_ratio = 1;

%% 4) Numerics
par.jacobian           = 'smallR';  % 'smallR' (as in the paper) or 'finiteR'
par.f_Hz               = logspace(2,8,801);  % analytic frequency grid [Hz]
par.sim.enabled        = true;      % stochastic simulation of the reduced model
par.sim.n_realizations = 10;        % averaged realizations (about 0.2 s each)
par.sim.dt             = 1e-3;      % time step, in units of t_s
par.sim.t_window       = [100,1000];% analysed time window, in units of t_s

%% Run the solver
res = reduced_model_noise_solver(par);

%% Plot
fig = plot_reduced_model_noise(res);

%% Save
here = fileparts(mfilename('fullpath'));
if isempty(here)
    here = pwd;
end
out_dir = fullfile(here,'output');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
name = strrep(sprintf('reduced_model_noise_P%.3g_delta%.3g_x%g', ...
    res.norm.P,res.norm.delta,par.noise_enhancement_ratio),'.','p');
exportgraphics(fig,fullfile(out_dir,[name,'.png']),'Resolution',300);
exportgraphics(fig,fullfile(out_dir,[name,'.pdf']),'ContentType','vector');
save(fullfile(out_dir,[name,'.mat']),'res');
fprintf('Saved %s.{png,pdf,mat} in %s\n',name,out_dir);
