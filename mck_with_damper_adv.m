function [x,a,diag,varargout] = mck_with_damper_adv(t,ag,M,C,K,k_sd,geom,orf,hyd,therm,num,cfg)
    nd = 1; if isfield(hyd,'n_parallel'), nd = max(1, round(hyd.n_parallel)); end
    Ao1 = max(orf.n_orf * (pi*geom.d_o^2/4), 1e-12);
    Ao  = nd * Ao1;              % toplam orifis alanı
    Ap1 = geom.Ap;
    Ap  = nd * Ap1;              % toplam piston alanı

    n = size(M,1); r = ones(n,1); Ns = n-1;
    z0 = zeros(2*n + 2*Ns + Ns + 2, 1);
    z0(2*n + (1:Ns)) = orf.p_amb;  z0(2*n + Ns + (1:Ns)) = orf.p_amb;
    z0(end-1)=therm.T0_C; z0(end)=therm.Ts0_C;

    % --- Ölçekli toleranslar (basınca ve akışa göre) ---
    rho0  = therm.rho_ref;
    dPcap = max(num.dP_cap, 1e5);
    Qcap_est = Ao * sqrt( 2*dPcap / max(rho0,100) );   % Ao = nd*Ao1

    AbsX = 1e-4;  AbsV = 1e-3;  AbsP = max(5e3, 1e-5 * dPcap);
    AbsQ = max(1e-5, 1e-3 * Qcap_est);  AbsT = 5e-3;
    AbsTol = [AbsX*ones(n,1); AbsV*ones(n,1); AbsP*ones(2*Ns,1); AbsQ*ones(Ns,1); AbsT*ones(2,1)];

    dt  = median(diff(t));
    opts = odeset('RelTol',7e-2,'AbsTol',AbsTol,'JPattern',local_JPattern(n), ...
                  'MaxStep',max(dt*15,3e-3),'InitialStep',max(dt*0.30,1.5e-3));

    agf = griddedInterpolant(t, ag, 'pchip', 'nearest');

    odef=@(tt,z) rhs(tt,z);
    sol=ode23tb(odef,[t(1) t(end)],z0,opts);

    t_end_sol=sol.x(end); t_use=t(t<=t_end_sol+1e-12); Z_use=deval(sol,t_use).';
    if t_end_sol < t(end)-1e-9
        opts2 = odeset(opts, 'RelTol',8e-2, 'MaxOrder', 2);
        try
            sol2=ode15s(odef,[t_end_sol t(end)],sol.y(:,end),opts2);
            t_use2=t(t>t_end_sol); Z_use2=deval(sol2,t_use2).'; Z=[Z_use;Z_use2];
            warning('ode23tb early stop at t=%.3f → continued to t=%.3f.',t_end_sol,t(end));
        catch
            Z=[Z_use; nan(numel(t)-numel(t_use), size(Z_use,2))];
            warning('Both solvers stopped early at t=%.3f s.',t_end_sol);
        end
    else
        Z=Z_use;
    end

    x = Z(:,1:n); v = Z(:,n+1:2*n);
    Fdev = dev_force_from_story(t, x, Z, k_sd, geom, hyd, therm, num, orf, cfg);
    a = ( -(M\(C*v.' + K*x.' + Fdev)).' - ag.*r.' );

    [drift, F_story, dP_orf_t, T_oil, T_steel, mu_t, E_cum, ...
     cav_frac_t, dP_q50, dP_q95, Q_abs_med, Q_abs_p95, ...
     cav_margin_t, cav_margin_min] = ...
        diagnostics_time_series(t, Z, n, Ns, k_sd, geom, orf, hyd, therm, num, cfg);

    diag = struct('drift',drift,'F_story',F_story,'dP_orf_time',dP_orf_t, ...
                  'dP_orf_time_max',max(dP_orf_t,[],2),'T_oil',T_oil,'T_steel',T_steel, ...
                  'mu_t',mu_t,'E_cum',E_cum,'cav_frac_t',cav_frac_t,'dP_q50',dP_q50,'dP_q95',dP_q95, ...
                  'Q_abs_med',Q_abs_med,'Q_abs_p95',Q_abs_p95, ...
                  'cav_margin_t',cav_margin_t,'cav_margin_min',cav_margin_min);

    if nargout>=4, varargout{1} = v; end

    % -------------------- iç fonksiyonlar --------------------
    function dz = rhs(tt,z)
        x  = z(1:n); v = z(n+1:2*n);
        p1 = z(2*n + (1:Ns)); p2 = z(2*n + Ns + (1:Ns));
        Q  = z(2*n + 2*Ns + (1:Ns)); T_o = z(end-1); T_s = z(end);

        mu_raw = therm.mu_ref * exp(therm.b_mu*(T_o - therm.T_ref_C));
        if cfg.use_thermal && cfg.on.mu_floor, mu=max(num.mu_min_phys,mu_raw); else, mu=cfg.use_thermal*mu_raw + (~cfg.use_thermal)*therm.mu_ref; end
        if cfg.use_thermal
            rho=max(100,therm.rho_ref/(1+therm.alpha_rho*(T_o-therm.T_ref_C)));
            beta=max(1e8,therm.beta0*exp(therm.b_beta*(T_o-therm.T_ref_C)));
            p_vap=p_vap_Antoine(T_o,therm,orf);
        else
            rho=therm.rho_ref; beta=therm.beta0; p_vap=p_vap_Antoine(therm.T_ref_C,therm,orf);
        end

        drift = x(2:end) - x(1:end-1);
        dvel  = v(2:end) - v(1:end-1);

        V1 = nd*hyd.V0 + Ap*drift;
        V2 = nd*hyd.V0 - Ap*drift;
        Vmin = nd * hyd.Vmin_fac * hyd.V0; V1=max(V1,Vmin); V2=max(V2,Vmin);

        % --- Orifis hidrolik kayıpları ---
        if cfg.use_orifice
            if cfg.on.CdRe
                Re = rho .* abs(Q) .* geom.d_o ./ max(Ao * mu, 1e-9);
                Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
                Cd = max(min(Cd, 1.2), 0.2);
            else
                Cd = orf.CdInf;
            end
            if cfg.on.Rkv, RQ = rho ./ max(2 * (Cd .* Ao).^2, 1e-12); else, RQ = 0 * Q; end

            Qcap = getfield_default(num,'Qcap_big', ...
                0.4 * ( max(max(orf.CdInf,orf.Cd0)*Ao, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
            if cfg.on.Qsat, Q_h = Qcap * tanh(Q ./ max(Qcap, 1e-9)); else, Q_h = Q; end

            dP_kv = RQ .* Q_h .* abs(Q_h);
            if cfg.on.Rlam
                R_lam  = (128 * mu * geom.Lori / (pi * geom.d_o^4)) / max(1, nd);
                dP_lam = R_lam .* Q;
            else
                dP_lam = 0 * Q;
            end
            dP_h = dP_lam + dP_kv;
        else
            dP_h = 0 * Q; Q_h = 0 * Q;
        end

        % --- Kavitasyon, dP_cap, atalet, leak ---
        if cfg.on.cavitation, p2_eff = max(p2, orf.cav_sf * p_vap); else, p2_eff = p2; end
        dP_raw = p1 - p2_eff;
        if cfg.on.dP_cap && isfinite(num.dP_cap)
            dP_eff = num.dP_cap * tanh(dP_raw ./ max(num.dP_cap, 1));
        else
            dP_eff = dP_raw;
        end
        if cfg.use_orifice && cfg.on.hyd_inertia, dQ = (dP_eff - dP_h) ./ max(hyd.Lh, 1e-12); else, dQ = 0 * Q; end
        if cfg.on.leak, Q_leak = hyd.K_leak * (p1 - p2); else, Q_leak = 0 * Q; end

        % --- Basınç ODE'leri ---
        if cfg.on.pressure_ode
            dp1 = (beta ./ V1) .* ( -Q - Q_leak - Ap * dvel );
            dp2 = (beta ./ V2) .* ( +Q + Q_leak + Ap * dvel );
        else
            dp1 = 0 * p1; dp2 = 0 * p2;
        end
        % Cavitation clamp
        if cfg.on.cavitation
            m1 = (p1 <= p_vap) & ((-Q - Q_leak - Ap*dvel) < 0);
            m2 = (p2 <= p_vap) & ((+Q + Q_leak + Ap*dvel) < 0);
            dp1(m1) = 0; dp2(m2) = 0;
        end

        % --- Damper kuvveti (yay + PF) — RESISTIVE-ONLY KLAMP BURADA ---
        w_pf = cfg.on.pressure_force * (1 - exp(-max(tt - cfg.PF.t_on, 0) / max(cfg.PF.tau, 1e-6)));
        F_story = k_sd * drift;

        if cfg.on.pressure_force
            dp_pf = (p1 - p2_eff);                % Ns x 1
            if isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
s = tanh(20*dvel);   % 20 ~ 1/(0.05 m/s) benzeri yumuşatma
                % İstersen daha yumuşak için: s = tanh(20*dvel);
                dp_pf = s .* max(0, s .* dp_pf);  % yalnız dirençli bileşen
            end
            F_story = F_story + Ap * (cfg.PF.gain * w_pf .* dp_pf);  % w_pf skaler
        end

        % --- Kuvvet dağıtımı, yapı ODE ---
        F        = zeros(n,1);
        F(1)     = -F_story(1);
        if n > 2, F(2:n-1) = F_story(1:end-1) - F_story(2:end); end
        F(n)     =  F_story(end);

        dv = M \ ( -C*v - K*x - F - M*r*agf(tt) );

        % --- Isı ODE'leri (toplam kayıp gücü) ---
        if cfg.use_thermal
            if cfg.use_orifice && cfg.on.Rlam, P_lam = dP_lam .* Q; else, P_lam = 0; end
            if cfg.use_orifice && cfg.on.Rkv,  P_kv  = dP_kv  .* Q; else, P_kv  = 0; end
            P_sum = sum(P_lam + P_kv); P_sum = max(P_sum, 0);

            dT_o = ( P_sum ...
                     - therm.hA_os   * (T_o - T_s) ...
                     - therm.hA_o_env* (T_o - therm.T_env_C) ) / max(therm.C_oil,   eps);
            dT_s = ( + therm.hA_os   * (T_o - T_s) ...
                     - therm.hA_s_env* (T_s - therm.T_env_C) ) / max(therm.C_steel, eps);
        else
            dT_o = 0; dT_s = 0;
        end

        dz = [ v; dv; dp1; dp2; dQ; dT_o; dT_s ];
    end

  function Jp = local_JPattern(nn)
    % Güvenli, tam dolu Jacobian şablonu (küçük sistemler için yeterli).
    % İsterseniz daha seyrek bir yapı kurabilirsiniz; bu sürüm garanti çalışır.
    Ns   = nn - 1;
    Ntot = 2*nn + 2*Ns + Ns + 2;  % 2n (x,v) + 2Ns (p1,p2) + Ns (Q) + 2 (T_o,T_s)
    Jp   = sparse(ones(Ntot, Ntot));
end

 end

function [drift, F_story, dP_orf_t, T_oil, T_steel, mu_t, E_cum, ...
          cav_frac_t, dP_q50, dP_q95, Q_abs_med, Q_abs_p95, ...
          cav_margin_t, cav_margin_min] = ...
          diagnostics_time_series(t, Z, n, Ns, k_sd, geom, orf, hyd, therm, num, cfg)

    nd  = 1; if isfield(hyd,'n_parallel'), nd = hyd.n_parallel; end
    Ao1 = orf.n_orf * (pi*geom.d_o^2/4);
    Ao  = nd * Ao1;
    Ap  = nd * geom.Ap;

    X  = Z(:,1:n);
    V  = Z(:,n+1:2*n);
    p1 = Z(:,2*n + (1:Ns));
    p2 = Z(:,2*n + Ns + (1:Ns));
    Q  = Z(:,2*n + 2*Ns + (1:Ns));
    T_oil   = Z(:,end-1);
    T_steel = Z(:,end);

    drift = X(:,2:end) - X(:,1:end-1);
    dvel  = V(:,2:end) - V(:,1:end-1);      % Nt x Ns

    mu_raw  = therm.mu_ref * exp(therm.b_mu*(T_oil - therm.T_ref_C));
    if cfg.use_thermal && cfg.on.mu_floor, mu_t = max(num.mu_min_phys, mu_raw);
    else,                                   mu_t = cfg.use_thermal.*mu_raw + (~cfg.use_thermal).*therm.mu_ref; end

    rho_t   = max(100, therm.rho_ref ./ (1 + therm.alpha_rho*(T_oil - therm.T_ref_C)));
    p_vap_t = p_vap_Antoine(T_oil, therm, orf);

    cav_margin_t  = min(p2 - p_vap_t, [], 2);
    cav_margin_min = min(cav_margin_t, [], 'omitnan');

    if cfg.on.cavitation, p2_eff = max(p2, orf.cav_sf * p_vap_t); else, p2_eff = p2; end

    % ---- PF kuvveti: RESISTIVE-ONLY burada uygulanıyor ----
    w_pf_vec = pf_weight(t, cfg) * cfg.PF.gain;
    if cfg.on.pressure_force
        dp_pf = (p1 - p2_eff);                         % Nt x Ns
        if isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
            s = tanh(20*dvel);   % 20 ~ 1/(0.05 m/s) benzeri yumuşatma
                            % Nt x Ns
            % İstersen: s = tanh(20*dvel);
            dp_pf = s .* max(0, s .* dp_pf);
        end
        F_story = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;  % Nt x Ns → (.*) yayımlı
    else
        F_story = k_sd*drift;
    end

    % ---- Hidrolik kayıplar (dP_h) + dP_orf_t ----
    if cfg.use_orifice
        if cfg.on.CdRe
            Re = rho_t .* abs(Q) .* geom.d_o ./ max(Ao .* mu_t, 1e-9);
            Cd = orf.CdInf - (orf.CdInf - orf.Cd0) ./ (1 + (Re./orf.Rec).^orf.p_exp);
            Cd = max(min(Cd, 1.2), 0.2);
        else
            Cd = orf.CdInf;
        end
        if cfg.on.Rkv, RQ = rho_t ./ max(2 * (Cd .* Ao).^2, 1e-12); else, RQ = 0 * Q; end

        Qcap = getfield_default(num,'Qcap_big', ...
            0.4 * ( max(max(orf.CdInf,orf.Cd0)*Ao, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
        if cfg.on.Qsat, Q_h = Qcap * tanh(Q ./ max(Qcap, 1e-9)); else, Q_h = Q; end

        dP_kv = RQ .* Q_h .* abs(Q_h);
        if cfg.on.Rlam
            R_lam_t = (128 * mu_t .* geom.Lori ./ (pi * geom.d_o^4)) / max(1, nd);
            dP_lam  = R_lam_t .* Q;
        else
            dP_lam  = 0 * Q;
        end
        dP_h = dP_lam + dP_kv;
    else
        Q_h = 0 * Q; dP_h = 0 * Q;
    end

    dP_raw = p1 - p2_eff;
    if cfg.on.dP_cap && isfinite(num.dP_cap)
        dP_eff = num.dP_cap * tanh( dP_raw ./ max(num.dP_cap, 1) );
    else
        dP_eff = dP_raw;
    end
    epsm     = max(1e3, double(num.softmin_eps));
    dP_orf_t = 0.5 * ( dP_eff + dP_h - sqrt( (dP_eff - dP_h).^2 + epsm^2 ) );

    cav_mask   = p2 < p_vap_t;
    cav_frac_t = mean(cav_mask, 2);

    if cfg.use_orifice
        P_lam = dP_lam .* Q; P_kv  = dP_kv  .* Q;
        P_sum = sum(P_lam + P_kv, 2); P_sum = max(P_sum, 0);
    else
        P_sum = zeros(size(t));
    end
    E_cum = cumtrapz(t, P_sum);

    dP_q50    = prctile(dP_orf_t, 50, 2);
    dP_q95    = prctile(dP_orf_t, 95, 2);
    qQ        = prctile(abs(Q(:)), [50 95]);
    Q_abs_med = qQ(1);
    Q_abs_p95 = qQ(2);
end

function F = dev_force_from_story(t, X, Z, k_sd, geom, hyd, therm, num, orf, cfg)
    Nt=size(X,1); n=size(X,2); Ns=n-1;
    drift = X(:,2:end) - X(:,1:end-1);
    V     = Z(:,n+1:2*n);
    dvel  = V(:,2:end) - V(:,1:end-1);    % Nt x Ns
    p1    = Z(:,2*n+(1:Ns)); 
    p2    = Z(:,2*n+Ns+(1:Ns)); 
    T_o   = Z(:,end-1);

    nd = 1; if isfield(hyd,'n_parallel'), nd = max(1, hyd.n_parallel); end
    Ap = nd * geom.Ap;

    if cfg.use_thermal, p_vap=p_vap_Antoine(T_o,therm,orf); else, p_vap=p_vap_Antoine(therm.T_ref_C,therm,orf); end
    if cfg.on.cavitation, p2_eff=max(p2,orf.cav_sf*p_vap); else, p2_eff=p2; end

    w_pf_vec = pf_weight(t,cfg) * cfg.PF.gain;

    if cfg.on.pressure_force
        dp_pf = (p1 - p2_eff);                  % Nt x Ns
        if isfield(cfg,'on') && isfield(cfg.on,'pf_resistive_only') && cfg.on.pf_resistive_only
s = tanh(20*dvel);   % 20 ~ 1/(0.05 m/s) benzeri yumuşatma
            % yumuşak istersen: s = tanh(20*dvel);
            dp_pf = s .* max(0, s .* dp_pf);
        end
        F_story = k_sd*drift + (w_pf_vec .* dp_pf) * Ap;
    else
        F_story = k_sd*drift;
    end

    F = zeros(n,Nt);
    F(1,:) = -F_story(:,1).';
    if n>2, F(2:n-1,:) = (F_story(:,1:end-1) - F_story(:,2:end)).'; end
    F(n,:) =  F_story(:,end).';
end
% ---- Buhar basıncı (Antoine) -----------------------------------------
function p_v = p_vap_Antoine(T_C, therm, ~)
    if isfield(therm,'antoine_A') && isfield(therm,'antoine_B') && isfield(therm,'antoine_C')
        A=therm.antoine_A; B=therm.antoine_B; C=therm.antoine_C;
    else, A=5.0; B=1700; C=-80; end
    T_C=double(T_C); p_v = 10.^(A - B./(C + T_C));
p_v = min(max(p_v, 5), 5e2);     % 5–500 Pa

end


