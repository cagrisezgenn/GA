function [f, aux] = eval_fitness_for_x(x, design_set, ...
    obj, cons, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, ...
    use_cache, fail_early_k)

    % ---- Önbellek anahtarı
    persistent memo
    if isempty(memo), memo = containers.Map('KeyType','char','ValueType','any'); end
    key = sprintf('set%d|%s', design_set, mat2str(x,8));

    if use_cache && isKey(memo,key)
        data = memo(key); f = data.f; aux = data.aux; return;
    end

    % ---- Tasarıma uygula ve değerlendir
    try
        % GA decode’u simulate içinde değil, doğrudan burada yapmaya gerek yok;
        % compute/evaluate fonksiyonları simulate’i çağırırken set/x’ı geçiriyoruz.

        % Amaç: J
        goal_src = tern(obj.use_scaled_for_goal,'scaled','raw');
        [J, ~] = compute_objective_over_records( ...
    goal_src, obj, tail_sec, ...
            t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
            t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
            M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
    design_set, x);

        % Kısıt: Penalty (erken çıkış desteği ile)
        cons_loc = cons; cons_loc.fail_early_k = fail_early_k;
        [Penalty, cons_detail] = evaluate_constraints_over_records( ...
            cons_loc, cons.src_for_constraints, obj, tail_sec, ...
            t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
            t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
            M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
    design_set, x);
      
% === J1 & J2 (split) hesap — Pareto günlüğü için ===
[J1_split, J2_split] = compute_objectives_split( ...
    tern(obj.use_scaled_for_goal,'scaled','raw'), obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num, ...
    cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
    design_set, x);

% aux içine koy
aux_J1 = J1_split; aux_J2 = J2_split;

% --- Pareto günlüğüne yaz ---
global PARETO;
PARETO.J1(end+1,1)  = aux_J1;
PARETO.J2(end+1,1)  = aux_J2;
PARETO.F(end+1,1)   = J + Penalty;
PARETO.Pen(end+1,1) = Penalty;
PARETO.set(end+1,1) = design_set;
PARETO.x{end+1,1}   = x(:).';
PARETO.feas(end+1,1)= (Penalty <= 1e-6);    % eşik: cezasız ≈ fizibıl

        f = J + Penalty;
        aux = struct('J',J,'Penalty',Penalty,'cons',cons_detail);

        if use_cache, memo(key) = struct('f',f,'aux',aux); end
    catch ME
    % Güvenli büyük ceza + ayrıntılı rapor
    f = 1e9;
    aux = struct('err', getReport(ME, 'extended', 'hyperlinks', 'off'));
end

end

