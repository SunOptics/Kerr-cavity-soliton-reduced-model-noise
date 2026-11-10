function fig = plot_reduced_model_noise(res,opts)
%PLOT_REDUCED_MODEL_NOISE Plot the output of reduced_model_noise_solver.
%
%   fig = plot_reduced_model_noise(res)
%   fig = plot_reduced_model_noise(res,opts)
%
%   Layout follows Fig. 4 of arXiv:2605.14614:
%     row (a), only if res.sim exists: one simulated trajectory of the peak
%       power |A|^2, Omega/(2*pi), Delta tau0/sqrt(Delta t) and phi/pi;
%     row (b): two-sided PSDs of the peak-intensity RIN, Omega/(2*pi),
%       tau0 and phi. Colored curves are the pathway contributions seeded
%       by the diffusion channels eta, Omega, tau0, phi; black is their
%       sum; cyan is the stochastic simulation.
%
%   opts (all optional):
%     visible     'on' (default) or 'off'
%     f_limits    frequency axis [Hz], default [1e2 1e8]
%     decades     y-span of the four PSD panels, default [8 7 14 7.5]
%     ylim        1x4 cell of explicit y-limits, overrides decades
%     sim_points  number of log-spaced points drawn of the simulated PSD
%     title       figure title; default lists the main parameters
%
%   Author:    Yifan Sun (yifan.sun@ulb.be)
%              Service OPERA-Photonics, Université libre de Bruxelles (ULB),
%              B-1050 Brussels, Belgium
%   Reference: Y. Sun et al., "Collective-Coordinate Fluctuations of
%              Driven-Dissipative Solitons", arXiv:2605.14614 (2026),
%              https://arxiv.org/abs/2605.14614, Fig. 4

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    defaults = struct('visible','on','f_limits',[1e2,1e8], ...
        'decades',[8,7,14,7.5],'sim_points',3000,'title',[]);
    defaults.ylim = {[],[],[],[]};
    names = fieldnames(opts);
    for kk = 1:numel(names)
        defaults.(names{kk}) = opts.(names{kk});
    end
    opts = defaults;

    ana = res.analytic;
    sim = res.sim;
    has_sim = ~isempty(sim);

    col = [0.85,0.33,0.40; ...   % eta
        0.15,0.50,0.85; ...      % Omega
        0.35,0.68,0.20; ...      % tau0
        0.98,0.52,0.18];         % phi
    col_sim = [0.36,0.72,0.66];

    nrow = 1+has_sim;
    fig = figure('Color','w','Visible',opts.visible, ...
        'Units','pixels','Position',[60,60,1500,320*nrow+70]);
    tl = tiledlayout(fig,nrow,4,'TileSpacing','compact', ...
        'Padding','compact');
    if isempty(opts.title)
        opts.title = sprintf(['L = %.4g m,  F = %.4g,  \\theta = %.3g,  ' ...
            'P_{in} = %.4g W  (P = %.4g),  \\delta = %.4g,  ' ...
            '\\tau_R = %.3g fs,  noise = %g \\times vacuum'], ...
            res.par.L_m,res.par.finesse,res.par.theta, ...
            res.norm.pump_power_W,res.norm.P,res.norm.delta, ...
            res.par.tauR_s*1e15,res.par.noise_enhancement_ratio);
    end
    title(tl,opts.title,'FontSize',11,'Interpreter','tex');

    if has_sim
        tr = sim.trajectory;
        t_ms = tr.t_s*1e3;
        ya = {tr.peak_power_W,tr.Omega_Hz*1e-9, ...
            tr.tau0_kick_sqrt_ps,tr.phi_rad/pi};
        la = {'$|\tilde A|^2$ (W)','$\tilde\Omega/(2\pi)$ (GHz)', ...
            '$\Delta\tilde\tau_0/\sqrt{\Delta\tilde t}$ ($\sqrt{\rm ps}$)', ...
            '$\tilde\phi/\pi$'};
        t_on_ms = sim.settings.t_noise_on*res.scales.t_s*1e3;
        for kk = 1:4
            ax = nexttile(tl);
            plot(ax,t_ms(2:end),ya{kk}(2:end),'Color',col_sim, ...
                'LineWidth',0.4);
            hold(ax,'on');
            xline(ax,t_on_ms,'k--','LineWidth',0.8);
            hold(ax,'off');
            xlim(ax,[0,t_ms(end)]);
            xlabel(ax,'$\tilde t$ (ms)','Interpreter','latex');
            ylabel(ax,la{kk},'Interpreter','latex');
            grid(ax,'on');
            panel_label_local(ax,sprintf('(a).%d',kk));
        end
    end

    ylab = {'RIN (1/Hz)', ...
        '$S_{\tilde\Omega\tilde\Omega}/(2\pi)^2$ (Hz$^2$/Hz)', ...
        '$S_{\tilde\tau_0\tilde\tau_0}$ (s$^2$/Hz)', ...
        '$S_{\tilde\phi\tilde\phi}$ (rad$^2$/Hz)'};
    f = ana.f_Hz;
    in_view = f >= opts.f_limits(1) & f <= opts.f_limits(2);
    handles = gobjects(1,6);
    for kk = 1:4
        ax = nexttile(tl);
        hold(ax,'on');
        if has_sim
            nf = numel(sim.f_Hz);
            idx = unique(round(logspace(0,log10(nf),opts.sim_points)));
            handles(1) = plot(ax,sim.f_Hz(idx),sim.psd(kk,idx), ...
                'Color',[col_sim,0.75],'LineWidth',2.2);
        end
        for cc = 1:4
            handles(2+cc) = plot(ax,f,squeeze(ana.contrib(kk,cc,:)), ...
                'Color',col(cc,:),'LineWidth',1.3);
        end
        handles(2) = plot(ax,f,ana.total(kk,:),'k','LineWidth',1.1);
        hold(ax,'off');
        set(ax,'XScale','log','YScale','log','Box','on');
        xlim(ax,opts.f_limits);
        if ~isempty(opts.ylim{kk})
            ylim(ax,opts.ylim{kk});
        else
            top = 10^ceil(log10(max(ana.total(kk,in_view)))+0.25);
            ylim(ax,[top/10^opts.decades(kk),top]);
        end
        xticks(ax,10.^(0:2:12));
        xlabel(ax,'$\tilde f$ (Hz)','Interpreter','latex');
        ylabel(ax,ylab{kk},'Interpreter','latex');
        grid(ax,'on');
        panel_label_local(ax,sprintf('(b).%d',kk));
    end
    labels = {'Simulation','Analytic total','\eta channel', ...
        '\Omega channel','\tau_0 channel','\phi channel'};
    keep = isgraphics(handles);
    lgd = legend(handles(keep),labels(keep),'FontSize',10);
    lgd.Layout.Tile = 'east';
end

function panel_label_local(ax,label)
    text(ax,0.97,0.93,label,'Units','normalized', ...
        'HorizontalAlignment','right','FontWeight','bold', ...
        'FontSize',11);
end
