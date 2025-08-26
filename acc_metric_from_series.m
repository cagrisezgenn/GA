function val = acc_metric_from_series(t, a, t5, t95, obj)
    % Arias penceresi içi metrikler
    if obj.use_arias_window
        w = (t>=t5 & t<=t95);
    else
        w = true(size(t));
    end
    aa = a(w); tt = t(w);
    if isempty(aa) || numel(aa)<2
        val = NaN; return;
    end
    switch lower(obj.acc_metric)
        case 'energy'
            val = trapz(tt, aa.^2);                 % enerji
        case 'rms+p95'
            rmsv = sqrt(mean(aa.^2,'omitnan'));
            p95  = prctile(abs(aa),95);
            val  = rmsv + obj.p95_penalty_w * p95;  % hibrit
        otherwise % 'rms'
            val = sqrt(mean(aa.^2,'omitnan'));
    end
end

