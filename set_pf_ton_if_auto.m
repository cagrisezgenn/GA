function cfg2 = set_pf_ton_if_auto(cfg2, t5, dt_on)
    %SET_PF_TON_IF_AUTO  Ensure cfg.PF.t_on if auto flag enabled.
    %   cfg2 = SET_PF_TON_IF_AUTO(cfg2, t5, dt_on) sets cfg2.PF.t_on to
    %   t5+dt_on when cfg2.PF.auto_t_on is true (default).  If auto_t_on is
    %   false, the existing cfg2.PF.t_on value is left untouched.  cfg2 is
    %   first passed through ensure_cfg_defaults to guarantee required
    %   fields.

    if nargin<3 || isempty(dt_on), dt_on = 0.5; end
    cfg2 = ensure_cfg_defaults(cfg2);
    if ~isfield(cfg2,'PF') || ~isfield(cfg2.PF,'auto_t_on') || cfg2.PF.auto_t_on
        cfg2.PF.t_on = t5 + dt_on;
    end
end
