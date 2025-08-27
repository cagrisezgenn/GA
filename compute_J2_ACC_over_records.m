function [J2, out] = compute_J2_ACC_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    if numel(varargin) >= 2, design_set = varargin{1}; x_ga = varargin{2};
    else, design_set = 0; x_ga = []; end
     
 
    useScaled = strcmpi(src,'scaled');

    useScaled = strcmpi(src,'scaled');
         useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R   = numel(tX);
    w_p = obj.p95_penalty_w;
    out = struct('a_rel', nan(R,1));

    mus = obj.mu_scenarios(:).';
    isWeighted = strcmpi(obj.mu_aggregate,'weighted');
    if isWeighted
        wmu = obj.mu_weights(:); wmu = wmu / max(sum(wmu), eps);
    end

    for r = 1:R
        % --- X yönü: A_env ref ve agg
        aref_X = NaN; aagg_X = NaN;
        if ~isempty(aX{r})
            [aref_X, aagg_X] = local_acc_ref_and_agg(tX{r}, aX{r}, t5x(r), t95x(r), 'X', r);
        end

        % --- Y yönü (varsa)
        aref_Y = NaN; aagg_Y = NaN;
        if ~isempty(aY{r})
            [aref_Y, aagg_Y] = local_acc_ref_and_agg(tY{r}, aY{r}, t5y(r), t95y(r), 'Y', r);
        end

        % --- Yön zarfı
        switch lower(obj.dir_mode)
            case 'xonly', a_ref = aref_X; a_agg = aagg_X;
            case 'yonly', a_ref = aref_Y; a_agg = aagg_Y;
            otherwise
                a_ref = max([aref_X, aref_Y], [], 'omitnan');
                a_agg = max([aagg_X, aagg_Y], [], 'omitnan');
        end

        out.a_rel(r) = a_agg / max(a_ref, eps);
    end

    J2 = cvar_from_samples(out.a_rel, obj.alpha_CVaR);

    % ---- local helper
    function [A_ref, A_agg] = local_acc_ref_and_agg(t, ag, t5, t95, dirStr, rid)
        % Referans (damper yok) mutlak ivme zarfı (RMS+p95)
        [~,a_rel0] = lin_MCK_consistent(t, ag, M, Cstr, K);  % relatif
        a_abs0 = a_rel0 + ag(:) * ones(1, size(a_rel0,2));  % her kat için mutlak
        w = (t >= t5 & t <= t95);
        A_i = zeros(1,size(a_abs0,2));
        for i=1:size(a_abs0,2)
            ai = a_abs0(w,i);
            rmsv = sqrt(mean(ai.^2,'omitnan'));
            p95  = prctile(abs(ai),95);
            A_i(i) = rmsv + w_p * p95;
        end
        A_ref = max(A_i);   % kat zarfı

        % Tasarım: kuyruk ekle + PF guard
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
                vals(k) = 5 * A_ref;  % fail durumunda güvenli büyük ceza
                continue;
            end
            a_relD = resp.y.a(1:numel(t), :);           % relatif, kuyruğu at
            a_absD = a_relD + ag(:) * ones(1,size(a_relD,2));

            Ai = zeros(1,size(a_absD,2));
            for i=1:size(a_absD,2)
                ai = a_absD(w,i);
                rmsv = sqrt(mean(ai.^2,'omitnan'));
                p95  = prctile(abs(ai),95);
                Ai(i)= rmsv + w_p * p95;
            end
            vals(k) = max(Ai);   % kat zarfı
        end

        if isWeighted, A_agg = nansum(wmu .* vals(:));
        else,          A_agg = max(vals, [], 'omitnan'); end
    end
end

