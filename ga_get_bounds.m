function [lb,ub,int_idx,names] = ga_get_bounds(set_id)
    switch set_id
        case 1
           %      d_o   Lori   mu_ref   Kd     d_w    D_m   n_turn    Dp    Lgap
            lb = [2e-3, 0.02, 1, 1.5e9, 3e-3, 0.06,   20, 0.10, 0.10];
            ub = [6e-3, 0.05, 1.50, 2.0e9,15e-3, 0.10,  40, 0.14, 0.24];
            int_idx = 7; % n_turn
            names = {'d_o','Lori','mu_ref','Kd','d_w','D_m','n_turn','Dp','Lgap'};

        case 2
            % [n_orf, Cd0, CdInf, Rec, p_exp, cav_sf, Lh, K_leak, resFactor]
            lb = [2, 0.55, 0.90, 2000, 1.0, 0.95, 1e-3, 1e-9, 10];
            ub = [2, 0.70, 1.10, 6000, 1.6, 1.00, 8e-3, 6e-9, 18];
            int_idx = 1; % n_orf
            names = {'n_orf','Cd0','CdInf','Rec','p_exp','cav_sf','Lh','K_leak','resFactor'};

        otherwise % 3  (n_orf SABİT=2)
            % [n_orf, Cd0, mu_ref, b_mu, beta0, b_beta, hA_os, dP_cap, Vmin_fac]
            lb = [2, 0.64, 0.6, -0.012, 1.2e9, -6e-3, 450, 3e8, 0.90,  10];
            ub = [2, 0.68, 1.2, -0.007, 2.2e9, -2e-3, 800, 6e8, 0.95, 18];
            int_idx = 1; % n_orf (sabit ama interface için sorun değil)
            names = {'n_orf','Cd0','mu_ref','b_mu','beta0','b_beta','hA_os','dP_cap','Vmin_fac'};

    end
end
