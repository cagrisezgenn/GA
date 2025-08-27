function cfg2 = cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg)
    % simulate() tasarımı içerde decode ediyor; burada yalnız cfg’yi aynen geçiriyoruz.
    %#ok<INUSD>
    cfg2 = ensure_cfg_defaults(cfg);
end
