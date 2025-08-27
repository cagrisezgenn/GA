function f = compact_log_wrapper(x, inner_fitfun)
% Tek satır log: [idx  #feval  J  (J+Penalty)  nViol]
% inner_fitfun: [f, aux] döndürebilir (f=J+Penalty), aux.J ve aux.cons.ratios içerebilir.
    persistent ROW FEVAL
    if isempty(ROW),   ROW   = 0; end
    if isempty(FEVAL), FEVAL = 0; end
    ROW   = ROW + 1;
    FEVAL = FEVAL + 1;

    J = NaN; nviol = 0;

    try
        [f_val, aux] = inner_fitfun(x);   % iki çıktı destekli
        if isstruct(aux)
            if isfield(aux,'J'), J = aux.J; end
            if isfield(aux,'cons') && isfield(aux.cons,'ratios') && isstruct(aux.cons.ratios)
                fn = fieldnames(aux.cons.ratios);
                vals = zeros(numel(fn),1);
                for i=1:numel(fn), vals(i) = aux.cons.ratios.(fn{i}); end
                nviol = sum(vals > 1+1e-12);
            end
        end
        f = f_val;

        % Hata raporu geldiyse konsola bas
        if isstruct(aux) && isfield(aux,'err') && ~isempty(aux.err)
            fprintf('ERR: %s\n', aux.err);
        end

    catch
        % inner_fitfun tek çıktı verirse
        try
            f = inner_fitfun(x);
        catch
            f = 1e9;   % güvenli büyük ceza
        end
    end

    fprintf('%6d %15d %13.3f %13.3f %8d\n', ROW, FEVAL, J, f, nviol);
end
