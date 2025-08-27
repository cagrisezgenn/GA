
function [Penalty, out] = evaluate_constraints_over_records( ...
    cons, src, obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    if numel(varargin) >= 2
        design_set = varargin{1}; x_ga = varargin{2};
    else
        design_set = 0; x_ga = [];
    end

      % Kaynak seçimi
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R   = numel(tX);
    mus = cons.mu_scenarios(:).';

    % ----- Erken çıkış sayaç ayarları -----
    fail_early_k = 0;
    if isfield(cons,'fail_early_k') && ~isempty(cons.fail_early_k) && isfinite(cons.fail_early_k)
        fail_early_k = max(0, round(cons.fail_early_k));
    end
    fail_count_global = 0;
    early_break_all   = false;

    % Toplayıcılar
    Fmax_records   = nan(R,1);
    stroke_records = nan(R,1);
    dpq_records    = nan(R,1);
    dT_records     = nan(R,1);
    cav_records    = nan(R,1);
    Qp95_records   = nan(R,1);
    any_fail_rec   = false(R,1);

    for r = 1:R
        % ---- Yön zarfı için X/Y ölçüleri
        rec_metrics = struct('Fmax',[],'stroke',[],'dpq',[],'dT',[],'cav95',[],'Qp95',[],'fail',[]);

        for dir = ["X","Y"]
            % Y yönü yoksa atla
            if dir=="Y" && (isempty(aY{r}) || all(isnan(aY{r})))
                continue;
            end

            % Kayıt/yön serileri
            if dir=="X"
                t = tX{r}; ag = aX{r}; t5=t5x(r); t95=t95x(r);
            else
                t = tY{r}; ag = aY{r}; t5=t5y(r); t95=t95y(r);
            end

            % ---- Tail ekle (simülasyon penceresi)
            dt    = median(diff(t));
            t_tail= (t(end)+dt : dt : t(end)+tail_sec).';
            t_s   = [t; t_tail];
            ag_s  = [ag; zeros(size(t_tail))];

            % ---- PF ramp koruması
cfg_dir = set_pf_ton_if_nan(cfg, t5, 0.5);


            % ---- μ-senaryosu zarfı
            Fmax_mu   = -Inf; stroke_mu = -Inf; dpq_mu = -Inf; dT_mu = -Inf; cav_mu = -Inf; Qp95_mu = -Inf;
            fail_mu   = false;

            % >>>>>>>>>>>>>>>>> μ DÖNGÜSÜ (BURADA) <<<<<<<<<<<<<<<<<
            for k = 1:numel(mus)
                resp = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                                M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir);
                                M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir, LOG);

                if ~resp.ok
                    fail_mu = true;
                    fail_count_global = fail_count_global + 1;

                    if fail_early_k > 0 && fail_count_global >= fail_early_k
                        early_break_all = true;   % tüm döngülerden çık
                        break;                    % μ-döngüsünü kır
                    end
                    continue;  % sıradaki μ
                end

                % ---- ölçümler (başarılı koşu)
                Fmax_mu   = max(Fmax_mu,   resp.F_max);
                stroke_mu = max(stroke_mu, resp.stroke_max);
                dpq_mu    = max(dpq_mu,    resp.dP_q_time(cons.dp.q));
                dT_mu     = max(dT_mu,     resp.dT_est);
                mask = (t_s >= t5) & (t_s <= t95);
    if any(mask)
        cav_p95_win = prctile(resp.diag.cav_frac_t(mask), 95);
    else
        cav_p95_win = prctile(resp.diag.cav_frac_t, 95);  % emniyetli geri dönüş
    end
    cav_mu = max(cav_mu, cav_p95_win);
                Qp95_mu   = max(Qp95_mu,   resp.metrics.Q_abs_p95);
            end
            
            % >>>>>>>>>>>>>>>>> μ DÖNGÜSÜ BİTTİ <<<<<<<<<<<<<<<<<<<<

            % μ-döngüsü eşik nedeniyle kırıldıysa, yön döngüsünü de bırak
            if early_break_all
                break;
            end

            % Bu yönün metriklerini biriktir
            rec_metrics.Fmax   = [rec_metrics.Fmax,   Fmax_mu];
            rec_metrics.stroke = [rec_metrics.stroke, stroke_mu];
            rec_metrics.dpq    = [rec_metrics.dpq,    dpq_mu];
            rec_metrics.dT     = [rec_metrics.dT,     dT_mu];
            rec_metrics.cav95  = [rec_metrics.cav95,  cav_mu];
            rec_metrics.Qp95   = [rec_metrics.Qp95,   Qp95_mu];
            rec_metrics.fail   = [rec_metrics.fail,   fail_mu];
        end

        % Yön döngüsünden erken çıktıysak kalan kayıtları INF kabul et
        if early_break_all
    any_fail_rec(r:end) = true;
    break;
end

if early_break_all
    Penalty = cons.pen.bigM * cons.pen.lambda.fail_bigM;
    out = struct();
    out.ratios = struct('spring_tau',0,'spring_slender',0,'stroke',0, ...
                        'force_cap',0,'dp_quant',0,'thermal_dT',0, ...
                        'cav_frac',0,'qsat_margin',0);
    out.any_fail = true;
    out.dpq_records    = dpq_records;
    out.Fmax_records   = Fmax_records;
    out.stroke_records = stroke_records;
    out.dT_records     = dT_records;
    out.cav_records    = cav_records;
    out.Qp95_records   = Qp95_records;
    return
end


        % ---- Yön zarfı (X,Y → max)
        Fmax_records(r)   = max(rec_metrics.Fmax,   [], 'omitnan');
        stroke_records(r) = max(rec_metrics.stroke, [], 'omitnan');
        dpq_records(r)    = max(rec_metrics.dpq,    [], 'omitnan');
        dT_records(r)     = max(rec_metrics.dT,     [], 'omitnan');
        cav_records(r)    = max(rec_metrics.cav95,  [], 'omitnan');
        Qp95_records(r)   = max(rec_metrics.Qp95,   [], 'omitnan');
        any_fail_rec(r)   = any(rec_metrics.fail);
    end


    % ---- Kayıt agregasyonları
    % Fmax ve stroke en-kötü kayıt
    Fmax_all   = max(Fmax_records,   [], 'omitnan');
    stroke_all = max(stroke_records, [], 'omitnan');
    dT_all     = max(dT_records,     [], 'omitnan');
    cav_all    = max(cav_records,    [], 'omitnan');
    Qp95_all   = max(Qp95_records,   [], 'omitnan');

    % Δp quantile: kayıtlar arası 'max' veya 'CVaR'
    switch lower(cons.dp.agg)
        case 'cvar'
            dpq_all = cvar_from_samples(dpq_records(:), cons.alpha_CVaR_cons);
        otherwise
            dpq_all = max(dpq_records, [], 'omitnan');
    end

    % ---- Oranlar (≥1 → ihlal)
    ratios = struct();

   % K1: yay kesme gerilmesi (yalnız yay kolu kuvveti)
if cons.on.spring_tau
    k_p_est = sh.G*sh.d_w^4 / (8*sh.n_turn*sh.D_m^3);           % [N/m]
    Fspring_max = k_p_est * stroke_all;                          % [N]
    ratios.spring_tau = spring_tau_ratio(Fspring_max, sh, cons.spring.tau_allow);
else
    ratios.spring_tau = 0;
end

    % K2: burkulma/serbest boy oranı
    if cons.on.spring_slender
        if strcmpi(cons.spring.L_free_mode,'fixed') && isfinite(cons.spring.L_free_fix)
            L_free = cons.spring.L_free_fix;
        else
            L_free = cons.spring.L_free_auto_fac * geom.Lgap; % yaklaşıklama
        end
        lambda = (L_free / max(sh.D_m,eps)) / max(cons.spring.lambda_max,eps);
        ratios.spring_slender = lambda;
    else
        ratios.spring_slender = 0;
    end

    % K3: strok
    if cons.on.stroke
        ratios.stroke = stroke_all / max(cons.stroke.util_factor*geom.Lgap, eps);
    else
        ratios.stroke = 0;
    end

    % K4: cihaz kuvvet limiti
    if cons.on.force_cap && isfinite(cons.force.F_cap)
        ratios.force_cap = Fmax_all / max(cons.force.F_cap, eps);
    else
        ratios.force_cap = 0;
    end

    % K5: Δp quantile
    if cons.on.dp_quant
        ratios.dp_quant = dpq_all / max(num.dP_cap, eps);
    else
        ratios.dp_quant = 0;
    end

    % K6: termal ΔT
    if cons.on.thermal_dT
        ratios.thermal_dT = dT_all / max(cons.thermal.cap_C, eps);
    else
        ratios.thermal_dT = 0;
    end

    % K7: kavitasyon
    if cons.on.cav_frac
        ratios.cav_frac = cav_all / max(cons.hyd.cav_frac_cap, eps);
    else
        ratios.cav_frac = 0;
    end

    % K8: Q satürasyon marjı
    if cons.on.qsat_margin
        Qcap = getfield_default(num,'Qcap_big', ...
    0.4 * ( max(max(orf.CdInf,orf.Cd0)*orf.Ao_eff, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
ratios.qsat_margin = Qp95_all / max(cons.hyd.Q_margin * Qcap, eps);

    else
        ratios.qsat_margin = 0;
    end

    % ---- Ceza hesabı
    hinge = @(r) max(0, r - 1);
    pwr   = cons.pen.power;

    pen = 0;
    if cons.on.spring_tau,     pen = pen + cons.pen.lambda.spring_tau     * hinge(ratios.spring_tau    )^pwr; end
    if cons.on.spring_slender, pen = pen + cons.pen.lambda.spring_slender * hinge(ratios.spring_slender)^pwr; end
    if cons.on.stroke,         pen = pen + cons.pen.lambda.stroke         * hinge(ratios.stroke        )^pwr; end
    if cons.on.force_cap && isfinite(cons.force.F_cap)
                               pen = pen + cons.pen.lambda.force_cap      * hinge(ratios.force_cap     )^pwr; end
    if cons.on.dp_quant,       pen = pen + cons.pen.lambda.dp_quant       * hinge(ratios.dp_quant      )^pwr; end
    if cons.on.thermal_dT,     pen = pen + cons.pen.lambda.thermal_dT     * hinge(ratios.thermal_dT    )^pwr; end
    if cons.on.cav_frac,       pen = pen + cons.pen.lambda.cav_frac       * hinge(ratios.cav_frac      )^pwr; end
    if cons.on.qsat_margin,    pen = pen + cons.pen.lambda.qsat_margin    * hinge(ratios.qsat_margin   )^pwr; end

    any_fail = any(any_fail_rec);
    if cons.on.fail_bigM && any_fail
        pen = pen + cons.pen.bigM * cons.pen.lambda.fail_bigM;
    end

    Penalty = pen;

    % Çıkış detayları
    out = struct();
    out.ratios    = ratios;
    out.any_fail  = any_fail;
    out.dpq_records = dpq_records;
    out.Fmax_records = Fmax_records;
    out.stroke_records = stroke_records;
    out.dT_records = dT_records;
    out.cav_records = cav_records;
    out.Qp95_records = Qp95_records;
end
