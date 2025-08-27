function v = cvar_from_samples(x, alpha)
    x = x(:); x = x(isfinite(x));
    if isempty(x), v = NaN; return; end
    q = quantile(x, 1-alpha);
    tail = x(x>=q);  % büyük-kötü kuyruk
    if isempty(tail), v = q; else, v = mean(tail); end
end
