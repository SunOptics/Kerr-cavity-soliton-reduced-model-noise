function res = reduced_model_noise_solver(par)
%REDUCED_MODEL_NOISE_SOLVER Noise spectra of a stationary cavity soliton
%from the four-coordinate reduced model x = (eta, Omega, tau0, phi).
%
%   res = reduced_model_noise_solver(par)
%
%   Model: Y. Sun et al., "Collective-Coordinate Fluctuations of
%   Driven-Dissipative Solitons", arXiv:2605.14614 (2026). Generalized
%   Lugiato-Lefever equation with Raman response, soliton ansatz
%       A = sqrt(2)*eta*sech(eta*(tau-tau0))*exp(-1i*Omega*(tau-tau0)+1i*phi).
%   Equation numbers in the comments refer to that paper.
%
%   Noise: white noise injected into the four collective coordinates,
%       D = r*D_vac,
%       D_vac = B*diag(eta/4, eta/12, pi^2/(48*eta^3), (1+pi^2/12)/(12*eta)),
%       B = 2*hbar*omega0/E_s                              [Eqs. (45), (51)]
%   r = par.noise_enhancement_ratio; r = 1 is the quantum-limited (vacuum)
%   noise. Every analytic PSD is proportional to r. Pump phase noise and
%   pump RIN are not included.
%
%   Input  par: struct of physical parameters; any missing field takes the
%          default listed in default_parameters_local below (the
%          fiber-ring example of Fig. 4 of the paper). See
%          run_reduced_model_noise.m for an annotated example.
%   Output res: struct with fields
%     scales      normalization scales and resonator quantities [SI]
%     norm        normalized parameters P, delta, d_R, B, d1
%     soliton     peak power, FWHM, energy, photon number, Raman shift
%     state       fixed point of the finite-R drift, Eq. (18)
%     M, D        Jacobian used for the analytic curves, diffusion matrix
%     eigenvalues eigenvalues of the finite-R Jacobian, stable (logical)
%     analytic    f_Hz, total (4 x Nf), contrib (4 x 4 x Nf: observable,
%                 seeding channel, frequency), two-sided PSDs in SI units
%     sim         stochastic reduced-model simulation (if enabled)
%   The four observables are: peak-intensity RIN [1/Hz],
%   S_OmOm/(2*pi)^2 [Hz^2/Hz], S_tau0tau0 [s^2/Hz], S_phiphi [rad^2/Hz].
%
%   Author:    Yifan Sun (yifan.sun@ulb.be)
%              Service OPERA-Photonics, Université libre de Bruxelles (ULB),
%              B-1050 Brussels, Belgium
%   Reference: Y. Sun, T. Bunel, S. Glazyrina, G. Semaan, F. Bretenaker,
%              S. Coen, S.-P. Gorza, and F. Leo, "Collective-Coordinate
%              Fluctuations of Driven-Dissipative Solitons",
%              arXiv:2605.14614 [physics.optics] (2026),
%              https://arxiv.org/abs/2605.14614
%   Figures:   the default parameters reproduce the reduced-model curves
%              of Fig. 4(a).1-4 and 4(b).1-4 (fiber ring, P = 3, delta = 9);
%              P = 4 with delta = 9 or 19 gives rows (c) and (d).
%              Model: Eqs. (5), (18), (23), (33)-(38), (51), (55)-(64).

    if nargin < 1 || isempty(par)
        par = struct();
    end
    par = merge_local(default_parameters_local(),par);
    validate_parameters_local(par);

    sc = physical_scales_local(par);
    nrm = normalized_parameters_local(par,sc);
    st = finiteR_state_local(nrm,sc,par);
    [M_finiteR,lambda,stable] = finiteR_jacobian_local(st);
    st_smallR = smallR_state_local(nrm.P,nrm.delta,nrm.dR);
    M_smallR = smallR_jacobian_local(st_smallR);
    nrm.d1 = st.d1;

    switch lower(par.jacobian)
        case 'smallr'
            ana_state = st_smallR;
            M = M_smallR;
        case 'finiter'
            ana_state = st;
            M = M_finiteR;
    end
    r = par.noise_enhancement_ratio;
    D = r*vacuum_diffusion_local(ana_state.eta,sc.B);
    ana = analytic_psd_local(M,D,ana_state,sc,par.f_Hz);

    res = struct();
    res.par = par;
    res.scales = sc;
    res.norm = nrm;
    res.soliton = soliton_properties_local(st,sc,lambda);
    res.state = st;
    res.state_smallR = st_smallR;
    res.M = M;
    res.M_finiteR = M_finiteR;
    res.M_smallR = M_smallR;
    res.eigenvalues = lambda;
    res.stable = stable;
    res.D = D;
    res.analytic = ana;
    res.sim = [];

    if par.verbose
        print_summary_local(res);
    end
    if ~stable
        warning('reducedModel:Unstable', ...
            ['The stationary soliton is linearly unstable (e.g. beyond the ' ...
            'Hopf boundary). The linear PSDs and the simulation are not ' ...
            'meaningful there.']);
    end

    if par.sim.enabled
        D_sim = r*vacuum_diffusion_local(st.eta,sc.B);
        res.sim = simulate_rm_local(st,D_sim,sc,par.sim);
        res.sim.variance_check = variance_check_local(M_finiteR,D_sim, ...
            res.sim,par.verbose);
    end
end

%% =====================================================================
%% Parameters
%% =====================================================================

function par = default_parameters_local()
% Fiber-ring example of Fig. 4 of the paper: P = 3, delta = 9.
    par = struct();
    par.L_m = 380;                     % roundtrip fiber length [m]
    par.beta1_s_per_m = 4.8823e-9;     % inverse group velocity [s/m]
    par.beta2_s2_per_m = -20e-27;      % group-velocity dispersion [s^2/m]
    par.gamma_per_W_m = 1.8e-3;        % Kerr coefficient [1/(W m)]
    par.finesse = 24;
    par.theta = 0.1;                   % power coupling ratio
    par.tauR_s = 3e-15;                % Raman time constant [s]
    par.lambda0_m = 1550e-9;           % pump wavelength [m]

    par.pump_power = 0.29512;           % see pump_power_unit
    par.pump_power_unit = 'W';         % 'W' or 'normalized' (P)
    par.detuning = 9;                  % see detuning_unit
    par.detuning_unit = 'normalized';  % 'normalized', 'rad' or 'Hz'

    par.noise_enhancement_ratio = 1;   % D = ratio*D_vac

    par.jacobian = 'smallR';           % 'smallR' (paper) or 'finiteR'
    par.f_Hz = logspace(2,8,801);      % frequencies of the analytic PSD
    par.verbose = true;

    par.sim.enabled = true;
    par.sim.dt = 1e-3;                 % normalized time step (units of t_s)
    par.sim.t_noise_on = 50;           % noise switched on at this time
    par.sim.t_window = [100,1000];     % analysis window (t1, t2]
    par.sim.n_realizations = 10;
    par.sim.seed = 1;
    par.sim.trajectory_t_max = 1000;   % stored span of realization 1
    par.sim.trajectory_stride = 50;    % keep every n-th step of that span
    par.sim.progress = true;
end

function validate_parameters_local(par)
    positive = {'L_m','beta1_s_per_m','gamma_per_W_m','finesse', ...
        'theta','lambda0_m','pump_power','noise_enhancement_ratio'};
    for kk = 1:numel(positive)
        value = par.(positive{kk});
        if ~isscalar(value) || ~isfinite(value) || value <= 0
            error('reducedModel:InvalidParameter', ...
                'par.%s must be a positive scalar.',positive{kk});
        end
    end
    if par.beta2_s2_per_m >= 0
        error('reducedModel:NormalDispersion', ...
            'Bright cavity solitons need anomalous dispersion, beta2 < 0.');
    end
    if par.tauR_s < 0
        error('reducedModel:InvalidParameter','par.tauR_s must be >= 0.');
    end
    if par.theta > 2*pi/par.finesse
        warning('reducedModel:Coupling', ...
            ['theta = %g exceeds the total roundtrip loss 2*pi/F = %.4g; ' ...
            'check the coupling ratio and the finesse.'], ...
            par.theta,2*pi/par.finesse);
    end
    if ~any(strcmpi(par.pump_power_unit,{'W','normalized'}))
        error('reducedModel:InvalidUnit', ...
            'par.pump_power_unit must be ''W'' or ''normalized''.');
    end
    if ~any(strcmpi(par.detuning_unit,{'normalized','rad','Hz'}))
        error('reducedModel:InvalidUnit', ...
            'par.detuning_unit must be ''normalized'', ''rad'' or ''Hz''.');
    end
    if ~any(strcmpi(par.jacobian,{'smallR','finiteR'}))
        error('reducedModel:InvalidJacobian', ...
            'par.jacobian must be ''smallR'' or ''finiteR''.');
    end
end

function sc = physical_scales_local(par)
    c0 = 299792458;
    hbar = 1.054571817e-34;

    sc = struct();
    sc.t_r = par.beta1_s_per_m*par.L_m;        % roundtrip time
    sc.FSR_Hz = 1/sc.t_r;
    sc.Omega_r = 2*pi*sc.FSR_Hz;
    sc.alpha_r = 2*pi/par.finesse;             % roundtrip energy loss, Eq. (7)
    sc.alpha = sc.alpha_r/sc.t_r;              % energy decay rate, Eq. (2)
    sc.cavity_linewidth_Hz = sc.alpha/(2*pi);  % = FSR/F
    sc.t_s = 2/sc.alpha;                       % Eq. (3)
    sc.f_s = 1/sc.t_s;
    sc.tau_s = sqrt(abs(par.beta2_s2_per_m)/(sc.alpha*par.beta1_s_per_m));
    sc.A_s2 = sc.alpha*par.beta1_s_per_m/(2*par.gamma_per_W_m);
    sc.A_s = sqrt(sc.A_s2);
    sc.E_s = sc.A_s2*sc.tau_s;
    sc.omega0 = 2*pi*c0/par.lambda0_m;
    sc.Q = sc.omega0/(2*sc.f_s);               % Eq. (7)
    sc.B = 2*hbar*sc.omega0/sc.E_s;            % Eq. (45)
    sc.hbar = hbar;
end

function nrm = normalized_parameters_local(par,sc)
% Eq. (6): P = 2*sqrt(theta*P_in)/(alpha*t_r*A_s), delta = 2*delta~/alpha_r.
    nrm = struct();
    switch lower(par.pump_power_unit)
        case 'w'
            nrm.pump_power_W = par.pump_power;
            nrm.P = 2*sqrt(par.theta*par.pump_power)/(sc.alpha_r*sc.A_s);
        case 'normalized'
            nrm.P = par.pump_power;
            nrm.pump_power_W = (nrm.P*sc.alpha_r*sc.A_s/2)^2/par.theta;
    end
    switch lower(par.detuning_unit)
        case 'normalized'
            nrm.delta = par.detuning;
        case 'rad'                             % roundtrip phase detuning
            nrm.delta = 2*par.detuning/sc.alpha_r;
        case 'hz'                              % detuning frequency
            nrm.delta = 2*pi*par.detuning*sc.t_s;
    end
    nrm.detuning_rad = nrm.delta*sc.alpha_r/2;
    nrm.detuning_Hz = nrm.delta/(2*pi*sc.t_s);
    nrm.dR = par.tauR_s/sc.tau_s;              % Eq. (6)
    nrm.B = sc.B;
    nrm.noise_enhancement_ratio = par.noise_enhancement_ratio;
end

%% =====================================================================
%% Reduced model: drift, stationary states, Jacobians, vacuum diffusion
%% =====================================================================

function G = drift_local(x,st)
% Finite-R drift, Eq. (18). Written without abs/conj so that it accepts a
% complex-step perturbation.
    eta = x(1);
    Om = x(2);
    phi = x(4);
    P = st.P;
    R = (pi/2)*Om/eta;                         % Eq. (19)
    S = sech(R);
    T = tanh(R);
    G = [ ...
        -2*eta+pi*P/sqrt(2)*S*cos(phi); ...
        -(16/15)*st.dR*eta^4-sqrt(2)*P*R*S*cos(phi); ...
        st.d1-2*Om+sqrt(2)*pi^2*P/(4*eta^2)*S*T*sin(phi); ...
        Om^2+eta^2-st.delta- ...
        3*pi*P/(4*sqrt(2)*eta)*S*sin(phi)*(1+(2/3)*R*T)];
end

function st = smallR_state_local(P,delta,dR)
% Small-R stationary relations, Eq. (23), on the upper branch phi >= 0
% (largest eta). delta decreases monotonically with phi on this branch.
    eta_of = @(phi) pi*P*cos(phi)/(2*sqrt(2));
    Om_of = @(phi) -(8/15)*dR*eta_of(phi).^4;
    mismatch = @(phi) eta_of(phi).^2+Om_of(phi).^2-1.5*tan(phi)-delta;
    phi = fzero(mismatch,[0,acos(1e-3)],optimset('TolX',1e-15));
    st = struct('P',P,'delta',delta,'dR',dR);
    st.eta = eta_of(phi);
    st.Omega = Om_of(phi);
    st.tau0 = 0;
    st.phi = phi;
    st.d1 = 2*st.Omega;                        % Eq. (23c)
    st.x = [st.eta;st.Omega;st.tau0;st.phi];
end

function st = finiteR_state_local(nrm,sc,par)
% Fixed point of Eq. (18): Newton on (G_eta, G_Omega, G_phi), which do not
% involve tau0 or d1, seeded by Eq. (23); then d1 from G_tau0 = 0.
    P = nrm.P;
    delta = nrm.delta;
    dR = nrm.dR;
    eta_max = pi*P/(2*sqrt(2));
    delta_max = eta_max^2+((8/15)*dR*eta_max^4)^2;
    if delta > delta_max
        P_min = 2*sqrt(2)*sqrt(delta)/pi;
        P_in_min = (P_min*sc.alpha_r*sc.A_s/2)^2/par.theta;
        error('reducedModel:NoSoliton', ...
            ['No stationary soliton: normalized delta = %.4g exceeds ' ...
            'delta_max = %.4g for P = %.4g (Eq. (23) at phi = 0).\n' ...
            'Reduce the detuning below %.4g (normalized) = %.4g Hz, or ' ...
            'raise the pump power above about %.4g W.'], ...
            delta,delta_max,P,delta_max,delta_max/(2*pi*sc.t_s),P_in_min);
    end
    seed = smallR_state_local(P,delta,dR);
    st = struct('P',P,'delta',delta,'dR',dR,'d1',0);
    z = [seed.eta;seed.Omega;seed.phi];
    h = 1e-20;
    for it = 1:50
        g = locked_drift_local(z,st);
        J = zeros(3);
        for jj = 1:3
            dz = zeros(3,1);
            dz(jj) = 1i*h;
            J(:,jj) = imag(locked_drift_local(z+dz,st))/h;
        end
        step = J\g;
        z = z-step;
        if norm(step) < 1e-15*max(1,norm(z))
            break;
        end
    end

    st.eta = z(1);
    st.Omega = z(2);
    st.tau0 = 0;
    st.phi = z(3);
    st.R = (pi/2)*st.Omega/st.eta;
    st.d1 = 2*st.Omega-sqrt(2)*pi^2*P/(4*st.eta^2)* ...
        sech(st.R)*tanh(st.R)*sin(st.phi);
    st.x = [st.eta;st.Omega;st.tau0;st.phi];
    st.residual = drift_local(st.x,st);
    if ~(st.eta > 0) || norm(st.residual) > 1e-10
        error('reducedModel:StationaryFailed', ...
            'Stationary solve failed: eta = %g, |G| = %g.', ...
            st.eta,norm(st.residual));
    end
end

function g = locked_drift_local(z,st)
    G = drift_local([z(1);z(2);0;z(3)],st);
    g = G([1,2,4]);
end

function [M,lambda,stable] = finiteR_jacobian_local(st)
% Complex-step Jacobian of Eq. (18), exact to round-off.
    h = 1e-20;
    M = zeros(4);
    for jj = 1:4
        dx = zeros(4,1);
        dx(jj) = 1i*h;
        M(:,jj) = imag(drift_local(st.x+dx,st))/h;
    end
    lambda = eig(M);
    [~,neutral] = min(abs(lambda));            % timing (translation) mode
    others = lambda;
    others(neutral) = [];
    stable = all(real(others) < 0);
end

function M = smallR_jacobian_local(st)
% Sparse small-R Jacobian, Eqs. (33)-(38).
    eta = st.eta;
    tp = tan(st.phi);
    m14 = -2*eta*tp;
    m21 = -(64/15)*st.dR*eta^3;
    m32 = pi^2*tp/(2*eta^2)-2;
    m41 = 2*eta+3*tp/(2*eta);
    m42 = -(16/15)*st.dR*eta^4;
    M = [-2,0,0,m14; ...
        m21,-2,0,0; ...
        0,m32,0,0; ...
        m41,m42,0,-3/2];
end

function D = vacuum_diffusion_local(eta,B)
% Eqs. (50)-(51), small-R vacuum diffusion matrix.
    D = B*diag([eta/4,eta/12,pi^2/(48*eta^3),(1+pi^2/12)/(12*eta)]);
end

function sol = soliton_properties_local(st,sc,lambda)
    sol = struct();
    sol.peak_power_W = 2*st.eta^2*sc.A_s2;           % Eq. (14)
    sol.FWHM_s = 2*acosh(sqrt(2))*sc.tau_s/st.eta;   % sech^2 intensity
    sol.energy_J = 4*st.eta*sc.E_s;                  % Eq. (12)
    sol.photon_number = sol.energy_J/(sc.hbar*sc.omega0);
    sol.raman_shift_Hz = st.Omega/(2*pi*sc.tau_s);
    sol.internal_mode_Hz = unique(abs(imag( ...
        lambda(abs(imag(lambda)) > 1e-9))))/(2*pi*sc.t_s);
end

%% =====================================================================
%% Analytic PSD and pathway decomposition
%% =====================================================================

function ana = analytic_psd_local(M,D,st,sc,f)
    f = f(:).';
    w = 2*pi*f*sc.t_s;                         % normalized angular frequency
    dk = diag(D).';
    contrib = zeros(4,4,numel(f));             % [observable, channel, f]
    for nn = 1:numel(f)
        H = (-1i*w(nn)*eye(4)-M)\eye(4);       % Eq. (55)
        contrib(:,:,nn) = abs(H).^2.*dk;       % |H_ik|^2 D_kk, Eq. (57)
    end

    % Normalized two-sided PSD -> physical units, Eqs. (59)-(60).
    scale = sc.t_s*[4/st.eta^2;1/(2*pi*sc.tau_s)^2;sc.tau_s^2;1];
    ana = struct();
    ana.f_Hz = f;
    ana.contrib = contrib.*scale;
    ana.total = reshape(sum(ana.contrib,2),4,[]);
    ana.observables = {'RIN','Omega/(2pi)','tau0','phi'};
    ana.units = {'1/Hz','Hz^2/Hz','s^2/Hz','rad^2/Hz'};
    ana.channels = {'eta','Omega','tau0','phi'};
    ana.rf_timing_phase_total = sc.Omega_r^2*ana.total(3,:);  % Eq. (64)
    ana.convention = 'two-sided, per Hz';
end

%% =====================================================================
%% Stochastic reduced-model simulation
%% =====================================================================

function sim = simulate_rm_local(st,D,sc,s)
% Euler-Maruyama of dx/dt = G(x) + F(t) with the finite-R drift, Eq. (18).
    dt = s.dt;
    n_on = round(s.t_noise_on/dt);
    n1 = round(s.t_window(1)/dt);
    n2 = round(s.t_window(2)/dt);
    if n1 < n_on
        warning('reducedModel:WindowBeforeNoise', ...
            'The analysis window starts before the noise is switched on.');
    end
    N = n2-n1;
    dts = dt*sc.t_s;                           % physical step [s]
    T = N*dts;
    nf = floor(N/2);
    f = (1:nf)/T;

    n_traj = min(round(s.trajectory_t_max/dt),n2);
    stride = max(1,round(s.trajectory_stride));
    n_tr = floor(n_traj/stride)+1;
    traj = zeros(n_tr,4);
    traj(1,:) = [st.eta,st.Omega,0,st.phi];

    sig = sqrt(diag(D)*dt);
    P = st.P;
    c_eta = pi*P/sqrt(2);
    c_R = (16/15)*st.dR;
    c_Om = sqrt(2)*P;
    c_tau = sqrt(2)*pi^2*P/4;
    c_phi = 3*pi*P/(4*sqrt(2));
    delta = st.delta;
    d1 = st.d1;
    kick_scale = sc.tau_s/sqrt(dts)/1e-6;      % Delta tau0/sqrt(Delta t) [sqrt(ps)]

    S = zeros(4,nf);
    var_sum = zeros(3,1);
    tic;
    for rr = 1:s.n_realizations
        rng(s.seed+rr-1,'twister');
        W = randn(4,n2-n_on);
        eta = st.eta;
        Om = st.Omega;
        tau = 0;
        phi = st.phi;
        Xe = zeros(N+1,1);
        Xo = zeros(N+1,1);
        Xt = zeros(N+1,1);
        Xp = zeros(N+1,1);
        keep_traj = rr == 1;
        for nn = 1:n2
            R = (pi/2)*Om/eta;
            Sr = 1/cosh(R);
            Tr = tanh(R);
            cp = cos(phi);
            sp = sin(phi);
            g1 = -2*eta+c_eta*Sr*cp;
            g2 = -c_R*eta^4-c_Om*R*Sr*cp;
            g3 = d1-2*Om+c_tau/eta^2*Sr*Tr*sp;
            g4 = Om^2+eta^2-delta-c_phi/eta*Sr*sp*(1+(2/3)*R*Tr);
            tau_old = tau;
            eta = eta+g1*dt;
            Om = Om+g2*dt;
            tau = tau+g3*dt;
            phi = phi+g4*dt;
            if nn > n_on
                mm = nn-n_on;
                eta = eta+sig(1)*W(1,mm);
                Om = Om+sig(2)*W(2,mm);
                tau = tau+sig(3)*W(3,mm);
                phi = phi+sig(4)*W(4,mm);
            end
            if nn >= n1
                kk = nn-n1+1;
                Xe(kk) = eta;
                Xo(kk) = Om;
                Xt(kk) = tau;
                Xp(kk) = phi;
            end
            if keep_traj && nn <= n_traj && mod(nn,stride) == 0
                traj(nn/stride+1,:) = [eta,Om, ...
                    (tau-tau_old)*kick_scale,phi];
            end
        end
        if any(~isfinite(Xe)) || any(Xe <= 0)
            ratio = 4*D(1,1)/(sc.B*st.eta);    % D(1,1) = ratio*B*eta/4
            error('reducedModel:SimulationDiverged', ...
                ['The stochastic simulation left the soliton branch; the ' ...
                'noise (noise_enhancement_ratio = %g) is too strong for ' ...
                'the reduced model.'],ratio);
        end

        I0 = 2*Xe(2:end).^2;                   % peak intensity / A_s^2
        y = [I0/mean(I0)-1, ...                % Delta I0 / I0
            Xo(2:end)/(2*pi*sc.tau_s), ...     % Omega~/(2 pi) [Hz]
            diff(Xt)*sc.tau_s/dts, ...         % kappa_n, Eq. (62)
            Xp(2:end)];                        % phi [rad]
        y = y-mean(y,1);
        Y = fft(y)*dts;
        S = S+(abs(Y(2:nf+1,:)).^2/T).';       % two-sided periodogram, Eq. (61)
        var_sum = var_sum+[var(Xe,1);var(Xo,1);var(Xp,1)];
        if s.progress
            fprintf('  simulation %d/%d done (%.1f s)\n', ...
                rr,s.n_realizations,toc);
        end
    end
    S = S/s.n_realizations;
    S(3,:) = dts^2*S(3,:)./(4*sin(pi*f*dts).^2);   % Eq. (63)

    sim = struct();
    sim.f_Hz = f;
    sim.psd = S;
    sim.observables = {'RIN','Omega/(2pi)','tau0','phi'};
    sim.variance_eta_Omega_phi = var_sum/s.n_realizations;
    sim.trajectory.t_s = (0:n_tr-1).'*stride*dts;  % physical time [s]
    sim.trajectory.peak_power_W = 2*traj(:,1).^2*sc.A_s2;
    sim.trajectory.Omega_Hz = traj(:,2)/(2*pi*sc.tau_s);
    sim.trajectory.tau0_kick_sqrt_ps = traj(:,3);
    sim.trajectory.phi_rad = traj(:,4);
    sim.settings = s;
    sim.T_s = T;
    sim.dt_s = dts;
end

function chk = variance_check_local(M,D,sim,verbose)
% Stationary covariance of the damped (eta, Omega, phi) subsystem of the
% simulated drift from A*Sigma + Sigma*A.' + D = 0 (column 3 of M is 0).
    idx = [1,2,4];
    A = M(idx,idx);
    Q = D(idx,idx);
    Sigma = reshape(-(kron(eye(3),A)+kron(A,eye(3)))\Q(:),3,3);
    chk = struct();
    chk.analytic = diag(Sigma);
    chk.simulated = sim.variance_eta_Omega_phi;
    chk.relative_error = chk.simulated./chk.analytic-1;
    if verbose
        fprintf(['  variance check (sim/analytic - 1): eta %+.3f, ' ...
            'Omega %+.3f, phi %+.3f\n'],chk.relative_error);
    end
end

%% =====================================================================
%% Utilities
%% =====================================================================

function print_summary_local(res)
    par = res.par;
    sc = res.scales;
    nrm = res.norm;
    st = res.state;
    sol = res.soliton;
    fprintf('\n=== Reduced-model soliton noise ===\n');
    fprintf(['Resonator: L = %.4g m, FSR = %.5g MHz, F = %.4g, ' ...
        'linewidth = %.4g kHz, Q = %.3g\n'],par.L_m,sc.FSR_Hz*1e-6, ...
        par.finesse,sc.cavity_linewidth_Hz*1e-3,sc.Q);
    fprintf(['Scales:    t_s = %.5g us, tau_s = %.5g ps, A_s^2 = %.5g W, ' ...
        'E_s = %.5g pJ, B = %.4g\n'],sc.t_s*1e6,sc.tau_s*1e12,sc.A_s2, ...
        sc.E_s*1e12,sc.B);
    fprintf(['Drive:     P_in = %.5g W, theta = %.4g  ->  P = %.5g;   ' ...
        'delta = %.5g (= %.4g rad = %.5g kHz);   d_R = %.4g\n'], ...
        nrm.pump_power_W,par.theta,nrm.P,nrm.delta,nrm.detuning_rad, ...
        nrm.detuning_Hz*1e-3,nrm.dR);
    fprintf('Noise:     D = %g x D_vac\n',par.noise_enhancement_ratio);
    fprintf(['Soliton:   eta = %.6g, Omega = %.4g, phi = %.6g;  ' ...
        'peak %.4g W, FWHM %.4g ps, %.3g photons, Raman shift %.4g GHz\n'], ...
        st.eta,st.Omega,st.phi,sol.peak_power_W,sol.FWHM_s*1e12, ...
        sol.photon_number,sol.raman_shift_Hz*1e-9);
    fprintf('Stability: eig(M) = %s\n',strjoin(arrayfun(@(z) ...
        sprintf('%.4f%+.4fi',real(z),imag(z)),res.eigenvalues.', ...
        'UniformOutput',false),', '));
    if ~isempty(sol.internal_mode_Hz)
        fprintf('           internal mode at %.4g kHz\n', ...
            sol.internal_mode_Hz*1e-3);
    end
    if res.stable
        fprintf('           linearly stable (one neutral timing mode)\n');
    end
    fprintf('Analytic curves use the %s Jacobian.\n',par.jacobian);
end

function out = merge_local(out,override)
    names = fieldnames(override);
    for kk = 1:numel(names)
        key = names{kk};
        if isfield(out,key) && isstruct(out.(key)) && ...
                isstruct(override.(key))
            out.(key) = merge_local(out.(key),override.(key));
        else
            out.(key) = override.(key);
        end
    end
end
