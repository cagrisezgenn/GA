function [J, out] = compute_objective_over_records( ...
src, obj, tail_sec, ...
t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    if numel(varargin) >= 2
        design_set = varargin{1}; x_ga = varargin{2};
    else
        design_set = 0; x_ga = [];
    end
    ...
    % local_design_ratios_one_dir içinde:
    % Simülasyon penceresi: aynı veriye kuyruk ekleyelim (PF ramp guard uyumu)
    % Kaynak seçici
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R = numel(tX);
    out(R) = struct('d_rel',NaN,'a_rel',NaN,'J_r',NaN);

    % Önce REFERANSLAR: (damper yok, μ=1.0) → her kayıt & (X,Y mevcutsa) yön için sabitlenir
    ref = struct('X',struct('d',nan(R,1),'a',nan(R,1)), ...
                 'Y',struct('d',nan(R,1),'a',nan(R,1)));

    for r=1:R
        % X yönü referans
        if ~isempty(aX{r})
            [dref, aref] = local_ref_metrics(tX{r}, aX{r}, t5x(r), t95x(r), obj, M,Cstr,K, n);
            ref.X.d(r)=dref; ref.X.a(r)=aref;
        end
        % Y yönü referans (varsa)
        if ~isempty(aY{r})
            [dref, aref] = local_ref_metrics(tY{r}, aY{r}, t5y(r), t95y(r), obj, M,Cstr,K, n);
            ref.Y.d(r)=dref; ref.Y.a(r)=aref;
        end
    end

    % Sonra TASARIM: her kayıt → yön zarfı → μ-agg → J_r
    for r=1:R
        [d_rel_X, a_rel_X] = local_design_ratios_one_dir('X', r);
        [d_rel_Y, a_rel_Y] = local_design_ratios_one_dir('Y', r);

        switch lower(obj.dir_mode)
            case 'xonly', d_rel = d_rel_X; a_rel = a_rel_X;
            case 'yonly', d_rel = d_rel_Y; a_rel = a_rel_Y;
            otherwise     % 'envelope'
                d_rel = max([d_rel_X, d_rel_Y], [], 'omitnan');
                a_rel = max([a_rel_X, a_rel_Y], [], 'omitnan');
        end

        % Ağırlıklı toplam
        J_r = obj.weights_da(1)*d_rel + obj.weights_da(2)*a_rel;

        out(r).d_rel = d_rel;
        out(r).a_rel = a_rel;
        out(r).J_r   = J_r;
    end

    % CVaR(α) hesap
    alpha = min(max(obj.alpha_CVaR,eps),0.99);
    Jlist = [out.J_r].';
    J = cvar_from_samples(Jlist, alpha);

    % ---- iç yardımcılar ----
 
    function [d_ref, a_ref] = local_ref_metrics(t, ag, t5, t95, obj, M,C,K, n)
        % Kuyruk eklemeden baz referans (damper yok)
        [x0,a0] = lin_MCK_consistent(t, ag, M, C, K);
        d_ref = max(abs(x0(:,min(obj.idx_disp_story,n))));
        a_ref = acc_metric_from_series(t, a0(:,min(obj.idx_acc_story,n)), t5, t95, obj);
        d_ref = max(d_ref, eps);
        a_ref = max(a_ref, eps);
    end

    function [d_rel_dir, a_rel_dir] = local_design_ratios_one_dir(which, r)
       
    switch upper(string(which))
        case "X"
            t=tX{r}; ag=aX{r}; t5=t5x(r); t95=t95x(r);
            d_ref=ref.X.d(r); a_ref=ref.X.a(r);
        case "Y"
            t=tY{r}; ag=aY{r}; t5=t5y(r); t95=t95y(r);
            d_ref=ref.Y.d(r); a_ref=ref.Y.a(r);
    end
    if isempty(ag) || ~isfinite(d_ref) || ~isfinite(a_ref)
        d_rel_dir = NaN; a_rel_dir = NaN; return;
    end


 
    % --- Simülasyon penceresi: aynı veriye kuyruk ekleyelim ---
 
    dt = median(diff(t));
tail_sec_loc = tail_sec;
    t_tail = (t(end)+dt : dt : t(end)+tail_sec_loc).';
    t_s    = [t;  t_tail];
    ag_s   = [ag; zeros(size(t_tail))];

    % PF ramp t_on: Arias t5 + 3
    cfg_dir = set_pf_ton_if_nan(cfg, t5, 0.5);


    % μ senaryoları
    mus   = obj.mu_scenarios(:).';
    d_vals = nan(size(mus));  a_vals = nan(size(mus));

    for k = 1:numel(mus)
        resp = simulate( ...
            design_set, x_ga, mus(k), t_s, ag_s, ...
            M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg_dir, LOG);

        if ~resp.ok
    fail_ratio = 5;              % istersen 3–10 arası dene
    d_vals(k) = fail_ratio * d_ref;
    a_vals(k) = fail_ratio * a_ref;
    continue;
end

        x  = resp.y.x(1:numel(t), :);   % kuyruğu at
        aR = resp.y.a(1:numel(t), :);

        d  = max(abs(x(:,min(obj.idx_disp_story,n))));
        am = acc_metric_from_series(t, aR(:,min(obj.idx_acc_story,n)), t5, t95, obj);

        d_vals(k) = d;  a_vals(k) = am;
    end

    % μ-aggregation
    switch lower(obj.mu_aggregate)
        case 'weighted'
            w = obj.mu_weights(:); w = w/sum(w);
            d_agg = nansum(w.*d_vals(:));
            a_agg = nansum(w.*a_vals(:));
        otherwise
            d_agg = max(d_vals, [], 'omitnan');
            a_agg = max(a_vals, [], 'omitnan');
    end

    d_rel_dir = d_agg / max(d_ref, eps);
    a_rel_dir = a_agg / max(a_ref, eps);
end

end

