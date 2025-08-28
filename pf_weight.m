function w = pf_weight(t, cfg)
%PF_WEIGHT Time-dependent pressure-force ramp factor.
%   w = PF_WEIGHT(t, cfg) returns the ramp weight for pressure-force
%   activation. It supports scalar or vector inputs t. The weight ramps
%   from 0 to 1 starting at cfg.PF.t_on with time constant cfg.PF.tau.
%
%   The result is multiplied by cfg.on.pressure_force so that the weight
%   is zero when the pressure-force mechanism is disabled.

    dt = max(t - cfg.PF.t_on, 0);
    w = cfg.on.pressure_force * (1 - exp(-dt ./ max(cfg.PF.tau, 1e-6)));
end
