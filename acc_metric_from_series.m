function val = acc_metric_from_series(t, a, t5, t95, obj)
    w = (t >= t5 & t <= t95);
    ai = a(w);
    switch lower(obj.acc_metric)
        case 'energy'
            val = sqrt(trapz(t(w), ai.^2));
        case 'rms+p95'
            rmsv = sqrt(mean(ai.^2,'omitnan'));
            p95  = prctile(abs(ai),95);
            val  = rmsv + obj.p95_penalty_w * p95;
        otherwise % 'rms'
            val = sqrt(mean(ai.^2,'omitnan'));
    end
end
