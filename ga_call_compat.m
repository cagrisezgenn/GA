function [xbest, fbest, output, pop, scores, exitflag] = ga_call_compat(fhandle, lb, ub, IntCon, opts)
% GA için sürüm-uyumlu çağrı (6/4/2 çıktı destekler)
    nvars = numel(lb);
    try
        % Yeni sürümler (6 çıktı)
        [xbest, fbest, exitflag, output, pop, scores] = ...
            ga(fhandle, nvars, [], [], [], [], lb, ub, [], IntCon, opts);
    catch
        try
            % Orta sürümler (4 çıktı)
            [xbest, fbest, exitflag, output] = ...
                ga(fhandle, nvars, [], [], [], [], lb, ub, [], IntCon, opts);
            pop = []; scores = [];
        catch
            % Eski sürümler (2 çıktı)
            [xbest, fbest] = ga(fhandle, nvars, [], [], [], [], lb, ub, [], IntCon, opts);
            exitflag = []; output = struct(); pop = []; scores = [];
        end
    end
end
