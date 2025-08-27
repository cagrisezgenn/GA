function [J1, out] = compute_J1_IDR_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    % Optional GA args
    if numel(varargin) >= 2, design_set = varargin{1}; x_ga = varargin{2};
    else, design_set = 0; x_ga = []; end

    % Kaynak seçimi
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R  = numel(tX);
    Ns = n - 1;
    if isscalar(h_story_m), h_story = repmat(h_story_m, Ns, 1);
    else,                   h_story = h_story_m(:); end

    out = struct('d_rel', nan(R,1));
    mus = obj.mu_scenarios(:).';
    isWeighted = strcmpi(obj.mu_aggregate,'weighted');
    if isWeighted
        wmu = obj.mu_weights(:); wmu = wmu / max(sum(wmu), eps);
    end

    for r = 1:R
        % --- X yönü referans IDR
        dref_X = NaN; dagg_X = NaN;
        if ~isempty(aX{r})
            [dref_X, dagg_X] = local_idr_ref_and_agg(tX{r}, aX{r}, t5x(r), t95x(r), 'X', r);
        end

        % --- Y yönü (varsa)
        dref_Y = NaN; dagg_Y = NaN;
        if ~isempty(aY{r})
            [dref_Y, dagg_Y] = local_idr_ref_and_agg(tY{r}, aY{r}, t5y(r), t95y(r), 'Y', r);
        end

        % --- Yön zarfı
        switch lower(obj.dir_mode)
            case 'xonly', d_ref = dref_X; d_agg = dagg_X;
            case 'yonly', d_ref = dref_Y; d_agg = dagg_Y;
            otherwise
                d_ref = max([dref_X, dref_Y], [], 'omitnan');
                d_agg = max([dagg_X, dagg_Y], [], 'omitnan');
        end

        out.d_rel(r) = d_agg / max(d_ref, eps);
    end

    J1 = cvar_from_samples(out.d_rel, obj.alpha_CVaR);

    % ---- local helper
    function [d_ref, d_agg] = local_idr_ref_and_agg(t, ag, t5, t95, dirStr, rid)
        % Referans (damper yok, μ=1.0), pencere içinde IDR zarfı
        [x0,~] = lin_MCK_consistent(t, ag, M, Cstr, K);
        w = (t >= t5 & t <= t95);
        drift0 = x0(:,2:end) - x0(:,1:end-1);                       % Nt x Ns
        idr0   = abs(drift0(w,:)) ./ (ones(sum(w),1) * h_story(:)'); % Nt_w x Ns
        d_ref  = max(idr0,[],'all');                                 % skaler

        % Tasarım: kuyruk ekle + PF ramp guard
        dt     = median(diff(t));
        t_tail = (t(end)+dt : dt : t(end)+tail_sec).';
        t_s    = [t; t_tail];
        ag_s   = [ag; zeros(size(t_tail))];

        cfg_dir = set_pf_ton_if_nan(cfg, t5, 0.5);

        vals = nan(size(mus));
        for k = 1:numel(mus)
            resp = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                            M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir);
                            M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir, LOG);
            if ~resp.ok
                vals(k) = 5 * d_ref;  % fail durumunda güvenli büyük ceza
                continue;
            end
            xD   = resp.y.x(1:numel(t), :);     % kuyruğu at
            driftD = xD(:,2:end) - xD(:,1:end-1);
            idrD   = abs(driftD(w,:)) ./ (ones(sum(w),1) * h_story(:)');
            vals(k)= max(idrD,[],'all');
        end

        if isWeighted, d_agg = nansum(wmu .* vals(:));
        else,          d_agg = max(vals, [], 'omitnan'); end
    end
end
end