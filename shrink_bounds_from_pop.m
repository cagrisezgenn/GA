function [lb2,ub2] = shrink_bounds_from_pop(pop, scores, lb, ub, keep_top, buf)
    if isempty(pop) || isempty(scores)
        lb2 = lb; ub2 = ub; return;
    end
    [~,ix] = sort(scores(:),'ascend');                 % en iyi küçük
    K = max(1, ceil(keep_top * size(pop,1)));
    P = pop(ix(1:K), :);                               % elitler
    p10 = prctile(P,10,1);
    p90 = prctile(P,90,1);
    span = max(p90 - p10, 1e-12);
    lb2 = max(lb, p10 - buf.*span);
    ub2 = min(ub, p90 + buf.*span);
end