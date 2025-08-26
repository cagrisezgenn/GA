function [geom, sh, orf, hyd, therm, num, ga] = decode_design_apply(ga, geom, sh, orf, hyd, therm, num)
    % --- set_id'i güvene al (string/char gelirse sayıya çevir) ---
    set_id = ga.design_set;
    if isstring(set_id) || ischar(set_id)
        set_id = str2double(set_id);
    end
    if isempty(set_id) || isnan(set_id)
        set_id = 1;
    end
    ga.design_set = double(set_id);
    % GA kapalıysa hiçbir şeyi değiştirme
    if ~isfield(ga,'enable') || ~ga.enable || isempty(ga.x)
        [lb,ub,ii,nn] = ga_get_bounds(ga.design_set);
        ga.lb=lb; ga.ub=ub; ga.int_idx=ii; ga.names=nn; ga.x_use=[];
        return;
    end

    [lb,ub,int_idx,names] = ga_get_bounds(ga.design_set);
    x = ga.x(:);
    if numel(x)~=numel(lb)
        error('ga.x uzunluğu set-%d için %d olmalı.', ga.design_set, numel(lb));
    end
    % Kenetleme + tamsayı
x = x(:);
lb = lb(:);
ub = ub(:);
x = max(lb, min(ub, x));
    if numel(int_idx)>=1
        idx = int_idx;
        x(idx) = round(x(idx));
        x = max(lb, min(ub, x));
    end

    switch ga.design_set
        case 1
            % [d_o, Lori, mu_ref, Kd, d_w, D_m, n_turn, Dp, Lgap]
            geom.d_o  = x(1);
            geom.Lori = x(2);
            therm.mu_ref = x(3);
            geom.Kd   = x(4);
            sh.d_w    = x(5);
            sh.D_m    = x(6);
            sh.n_turn = x(7);
            geom.Dp   = x(8);
            geom.Lgap = x(9);
        case 2
            % [n_orf, Cd0, CdInf, Rec, p_exp, cav_sf, Lh, K_leak, resFactor]
            orf.n_orf   = x(1);
            orf.Cd0     = x(2);
            orf.CdInf   = max(x(3), x(2)+0.05); % tutarlılık
            orf.Rec     = x(4);
            orf.p_exp   = x(5);
            orf.cav_sf  = x(6);
            hyd.Lh      = x(7);
            hyd.K_leak  = x(8);
            therm.resFactor = x(9);
        otherwise % 3
            % [n_orf, Cd0, mu_ref, b_mu, beta0, b_beta, hA_os, dP_cap, Vmin_fac]
            orf.n_orf   = x(1);
            orf.Cd0     = x(2);
            therm.mu_ref= x(3);
            therm.b_mu  = x(4);
            therm.beta0 = x(5);
            therm.b_beta= x(6);
            therm.hA_os = x(7);
            num.dP_cap  = x(8);
            hyd.Vmin_fac= x(9);
therm.resFactor = x(10);
    end

    ga.lb=lb; ga.ub=ub; ga.int_idx=int_idx; ga.names=names; ga.x_use=x;
    % ... (x hazırlandıktan SONRA)
    xdisp = strjoin( compose('%.3g', x(:).'), ', ' );  % 9 sayı -> "a, b, c..."
 % (Sessiz mod) — yalnızca LOG.verbose_decode=true ise yaz
try
    if evalin('base','exist(''LOG'',''var'')') && evalin('base','isstruct(LOG)') ...
            && evalin('base','isfield(LOG,''verbose_decode'') && LOG.verbose_decode')
        fprintf('GA decode: set=%s | x_use = [%s]\n', num2str(ga.design_set), xdisp);
    end
catch
    % base workspace erişimi yoksa sessizce geç
end
end
