function cfg = ensure_cfg_defaults(cfg)
    if nargin==0 || ~isstruct(cfg), cfg = struct(); end
    % top-level toggles
    def_top = {'use_orifice',true; 'use_thermal',true};
    for i=1:size(def_top,1)
        if ~isfield(cfg,def_top{i,1}), cfg.(def_top{i,1}) = def_top{i,2}; end
    end
    % nested on-struct
    if ~isfield(cfg,'on') || ~isstruct(cfg.on), cfg.on = struct(); end
    def_on = {'CdRe',true; 'Rlam',true; 'Rkv',true; 'Qsat',true; 'cavitation',true; ...
              'dP_cap',true; 'hyd_inertia',true; 'leak',true; ...
              'pressure_ode',true; 'pressure_force',true; 'mu_floor',true; ...
              'pf_resistive_only', false};   % <<< YENİ: resistive-only anahtarı
    for i=1:size(def_on,1)
        f = def_on{i,1}; if ~isfield(cfg.on,f), cfg.on.(f) = def_on{i,2}; end
    end
    % PF struct
    if ~isfield(cfg,'PF') || ~isstruct(cfg.PF), cfg.PF = struct(); end
    if ~isfield(cfg.PF,'mode'), cfg.PF.mode = 'ramp'; end
    if ~isfield(cfg.PF,'t_on'), cfg.PF.t_on = NaN; end
    if ~isfield(cfg.PF,'tau'),  cfg.PF.tau  = 2.5; end
    if ~isfield(cfg.PF,'gain'), cfg.PF.gain = 0.6; end
end
