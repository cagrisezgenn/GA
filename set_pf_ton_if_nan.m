function cfg2 = set_pf_ton_if_nan(cfg2, t5, dt_on)
    if nargin<3 || isempty(dt_on), dt_on = 0.5; end
    cfg2 = ensure_cfg_defaults(cfg2);
    if ~isfield(cfg2,'PF') || ~isfield(cfg2.PF,'t_on') || isnan(cfg2.PF.t_on)
        cfg2.PF.t_on = t5 + dt_on;
    end
end
