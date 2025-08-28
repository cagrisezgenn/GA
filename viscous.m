%% =====================================================================
%  10-Katlı Çerçeve — ODE-only viskoz damper modeli (CLEAN / AHENK + Switch)
%  (Sıkışabilirlik + Hat ataletı + Cd(Re) orifis + kavitasyon + 2-düğümlü ısı)
%  Laminer kayıp hidrolikte: Δp_lam(T) = R_lam(T)*Q
%  Kavitasyon: p2_eff = max(p2, cav_sf * p_vap(T)) hem akışta hem kuvvette
%  Solver: ode23tb (+ güvenli deval + ode15s fallback)
%  NOT: Antoine katsayılarını yağınıza göre kalibre edin.
% =====================================================================

clear; clc; close all;

% Global log anahtarları
global LOG; LOG = struct('verbose_decode', false);
global PARETO;
global PSA_PREP_CACHE;
if isempty(whos('global','PARETO')) || isempty(PARETO)
    PARETO = struct('J1',[],'J2',[],'F',[],'Pen',[], ...
                    'set',[],'x',{{}},'feas',[]);
end

%% -------------------- Kullanıcı anahtarları ---------------------------

% (1) Kayıt/yön seçimi + görselleştirme
vis.make_plots          = true;
vis.sel_rec             = 1;         % 1..R
vis.sel_dir             = 'x';       % 'X' veya 'Y'
vis.dual_run_if_needed  = true;      % plot_source ≠ sim_source ise ikinci koşu yap

% (2) Örnekleme / yeniden örnekleme anahtarları
prep.target_dt      = 0.005;     % hedef dt (s)
prep.resample_mode  = 'auto';     % 'auto'|'off'|'force'
prep.regrid_method  = 'pchip';   % 'pchip'|'linear'
prep.tol_rel        = 1e-6;      % dt eşitlik toleransı (göreli)

% (3) Şiddet eşitleme / PSA / CMS anahtarları
pp.on.intensity      = true;       % band-ortalama Sa(+CMS) normalizasyonu (yalnız GA için)
pp.on.CMS            = false;      % cms_target.mat varsa ve kullanmak istersen true
pp.gammaCMS          = 0.50;       % hibrit ağırlık (0→yalnız band, 1→yalnız CMS)

pp.PSA.zeta          = 0.05;       % SDOF sönüm oranı
pp.PSA.band_fac      = [0.8 1.2];  % T1 bandı (T1±%20)
pp.PSA.Np_band       = 15;         % band içi periyot sayısı
pp.PSA.downsample_enabled = true;  % SA hesabı için downsample kullan (false → ham dt)
pp.PSA.downsample_dt      = 0.02;  % [s] downsample hedef dt; yalnız downsample_enabled=true iken
pp.PSA.use_parfor        = false;  % Parallel Toolbox varsa denersin

% (4) Arias penceresi & kuyruk
pp.on.arias          = true;
pp.tail_sec          = 20;

% (5) Simülasyon ve grafik için veri kaynağı seçimi
%     'raw'    : ham ivmeler (genlik korunur)
%     'scaled' : şiddet eşitlenmiş ivmeler (GA/normalize amaçlı)
sel.sim_source  = 'scaled';
sel.plot_source = 'raw';

%% -------------------- Amaç fonksiyonu anahtarları ---------------------

% (A) Hedef katlar ve metrik
obj.idx_disp_story   = 10;        % d_10 → tepe yer değiştirme
obj.idx_acc_story    = 3;         % a_3  → ivme metriği
obj.acc_metric       = 'rms+p95';     % 'rms' | 'energy' | 'rms+p95' (hibrit)
obj.p95_penalty_w    = 0.20;      % hibritte pik cezası ağırlığı

%% -------------------- Kısıt anahtarları -------------------------------

% Kısıtları aç/kapa ve eşik/ceza ayarları
cons.on.spring_tau     = true;    % K1: yay kesme gerilmesi (τ_max ≤ τ_allow)
cons.on.spring_slender = true;    % K2: L_free/D_m ≤ λ_max
cons.on.stroke         = true;    % K3: max|drift| ≤ 0.9*L_gap
cons.on.force_cap      = true;    % K4: max|F_story| ≤ F_cap
cons.on.dp_quant       = true;    % K5: q≈0.99 Δp_orf ≤ dP_cap
cons.on.thermal_dT     = true;    % K6: ΔT_est ≤ ΔT_cap
cons.on.cav_frac       = false;    % K7: kavitasyon payı sınırı
cons.on.qsat_margin    = true;    % K8: Q 95p ≤ margin*Qcap_big
cons.on.fail_bigM      = true;    % K9: Simülasyon başarısızsa büyük ceza

% Eşikler / limitler
cons.spring.tau_allow   = 300e6;  % [Pa] yay çeliği için tipik, gerekirse güncelle
cons.spring.lambda_max  = 12.0;   % boy/çap sınırı
% true -> serbest yay boyu sabit değerden alınır; false -> boşluk oranına göre hesaplanır
cons.spring.use_fixed_length = false;   % varsayılan: false
cons.spring.L_free_fixed     = NaN;     % [m] use_fixed_length=true ise kullanılan serbest boy
cons.spring.L_free_auto_fac = 2.2; % 'auto' modda L_free ≈ fac * L_gap

cons.stroke.util_factor = 0.90;   % izinli strok = 0.90*L_gap

cons.force.F_cap        = 2e6;    % [N] cihaz kuvvet sınırı (varsayılan: 2e6); kısıtı cons.on.force_cap ile aç/kapa
cons.dp.q               = 0.99;   % Δp_orf zaman-içi quantile
cons.dp.agg             = 'max';  % kayıtlar arası: 'max' | 'cvar'
cons.alpha_CVaR_cons    = 0.20;   % yalnız 'cvar' seçilirse kullanılır

cons.thermal.cap_C      = 30.0;   % [°C] yağ ısınma sınırı
cons.hyd.cav_frac_cap   = 0.15;   % kavitasyon zaman oranı sınırı (95p)
cons.hyd.Q_margin       = 0.90;   % Q 95p ≤ margin*Qcap_big

% Ceza ayarları
cons.pen.power  = 2;              % hinge^power
cons.pen.bigM   = 1e6;            % simülasyon başarısızlığında eklenecek ceza
cons.pen.lambda = struct( ...     % kısıt başına ağırlıklar
    'spring_tau',     0.3, ...%0.3
    'spring_slender', 0.2, ...%0.2
    'stroke',         1.2, ...%1.2
    'force_cap',      0, ...%0
    'dp_quant',       0.5, ...%0.5
    'thermal_dT',     0.2, ...%0.2
    'cav_frac',       1.5, ...%1.5
    'qsat_margin',    0.3, ...%0.3
    'fail_bigM',      1 );%

% Kısıt simülasyon kaynağı ve μ seti (amaçla tutarlı kalsın)
cons.src_for_constraints = 'raw';           % kısıtları 'raw' ile değerlendir
if isfield(obj,'mu_scenarios')
    cons.mu_scenarios = obj.mu_scenarios;
else
    cons.mu_scenarios = [0.75 1.00 1.25];  % varsayılan; aşağıda obj tanımlanınca senkronlayacağız
end

% (B) Arias penceresi içi ölçüm
obj.use_arias_window = true;      % true → [t5,t95] aralığında metrik
obj.window_source    = 'same';    % 'same' → kullanılan seri üzerinden hesapla

% (C) Yön zarfı
obj.dir_mode         = 'Xonly'; % 'envelope' (=max(X,Y)) | 'Xonly' | 'Yonly'

% (D) μ-robustluk
obj.mu_scenarios = [0.9 1.0 1.1];
obj.mu_aggregate = 'weighted';
obj.mu_weights   = [0.25 0.50 0.25];
% Kısıt μ-senaryolarını amaçla hizala
cons.mu_scenarios = obj.mu_scenarios;

% (E) Risk agregasyonu (kayıtlar üzerinde)
obj.use_scaled_for_goal = true;   % amaç için 'scaled' setiyle çalış (önerilir)
obj.alpha_CVaR       = 0.20;      % tail payı
obj.weights_da       = [0.5 0.5]; % [w_disp, w_acc] toplam 1.0

% (F) Referans tanımı
% d_ref ve a_ref: aynı kayıt/yön için (damper yok, μ=1.0), sabit baz koşusundan.
% NOT: Bu referanslar bir kez hesaplanır ve tasarımla değişmez.

%% -------------------- Model anahtarları (blok-bazlı aç/kapa) ----------

% NOT: Bu blok dosyanın BAŞINDA olmalı. Aşağıdaki varsayılanları
% dilediğin gibi değiştir; ayrıca ensure_cfg_defaults() eksikleri tamamlar.
cfg.use_orifice = true;
cfg.use_thermal = true;
cfg.on.CdRe            = true;
cfg.on.Rlam            = true;
cfg.on.Rkv             = true;
cfg.on.Qsat            = true;
cfg.on.cavitation      = true;
cfg.on.dP_cap          = true;
cfg.on.hyd_inertia     = true;
cfg.on.leak            = true;
cfg.on.pressure_ode    = true;
cfg.on.pressure_force  = true;
cfg.on.mu_floor        = true;


% Basınç-kuvvet rampa/kazanç (PF)
cfg.PF.mode       = 'ramp';
cfg.PF.auto_t_on  = true;   % true -> t_on = t5 + 0.5 via set_pf_ton_if_auto
cfg.PF.t_on       = 0;      % auto_t_on=false ise manuel değer (s)
cfg.PF.tau        = 0.9;    % 4.0 → 1.5  (önerim)
cfg.PF.gain       = 1.7;    % 0.45 → 1.0 (önerim)
cfg.on.pf_resistive_only = true;  % sadece rezistif PF bileşeni kullan

% Eksik/yanlış alanlar için güvenli tamamlayıcı (guard; PF.auto_t_on dahil)
cfg = ensure_cfg_defaults(cfg);

%% -------------------- GA/Opt tasarım anahtarları ----------------------

ga.enable      = false;    % GA tasarım vektörü uygula? (false → baz set)
ga.design_set  = 1;        % 1|2|3 (aşağıdaki set tanımları)
ga.x           = [];       % Örn: set-1 için 9x1 vektör, boşsa uygulanmaz
% Örnek: otomatik orta-nokta denemesi
% [lb,ub,~,~] = ga_get_bounds(ga.design_set); ga.x = 0.5*(lb+ub); ga.enable=true;

%% -------------------- GA ÇALIŞTIRICI (opsiyonel) ----------------------

% --- Değişken sınırları ve IntCon hazırlığı ---
[lb,ub,int_idx,~] = ga_get_bounds(ga.design_set);
if isempty(int_idx)
    IntCon = [];           % tamsayı değişken yok
else
    IntCon = int_idx(:)';  % tamsayı indeks vektörü
end

% GA anahtarları
ga.opt.enable       = true;
ga.opt.seed         = 11;
ga.opt.popsize      = 6;
% Çok aşamalı akış: önce set-1, sonra set-3
ga.opt.multi_stage  = [1];
% Aşama başına nesil sayısı (opsiyonel):
ga.opt.maxgen_stage = [2];  % her blok için nesil
ga.opt.maxgen       = ga.opt.maxgen_stage(1);   % ilk aşama için başlangıç
ga.opt.use_parallel  = true;      % GA.UseParallel
ga.opt.lhs_refine    = true;      % 2-aşamalı LHS (coarse→refine)
ga.opt.cache_enable  = true;      % aynı x için memoize
ga.opt.fail_early_k  = 6;         % erken çıkış eşiği (başarısız sim sayısı)
ga.opt.save_log      = 'runLog_ga.mat';  % nüfus & en iyi çözüm logu
% daraltma parametreleri
ga.opt.keep_top     = 0.20;   % en iyi %20'yi tut
ga.opt.buffer_fac   = 0.10;   % p10–p90 etrafına %10 tampon

%% -------------------- Girdiler (7 kayıt; sütunlar: t, ax, (ops) ay) ---

Sall = load('acc_matrix.mat');
fn   = fieldnames(Sall);
fn   = fn(startsWith(fn,'acc_matrix'));
R    = numel(fn);
if R==0, error('acc_matrix.mat içinde acc_matrix* isimli dizi bulunamadı.'); end

%% -------------------- Yapı (T1 için gerekli) --------------------------

n  = 10;
m  = 2.2e6 * ones(n,1);
k  = 2.95e8 * ones(n,1);
c0 = 2.55e6 * ones(n,1);
[M,K,Cstr] = make_KCM(n,m,k,c0);
[~,D] = eig(K,M); w = sqrt(sort(diag(D),'ascend')); T1 = 2*pi/w(1);
% --- IDR için kat yükseklik(leri) ---
h_story_m = 3.0 * ones(n-1,1);   % tüm katlar 3.0 m ise
% h_story_m = [h1; h2; ...; h_{n-1}];   % kat kat farklı ise (alternatif)


%% -------------------- PSA fonksiyon seçimi ----------------------------

f_band = @sdof_PSA_band_avg_aug;   % tek ODE, çok periyot
f_vec  = @sdof_PSA_vec_aug_ode;

%% -------------------- Kayıtları oku → RAW (tekil & hedef dt kontrolü) -

t_rawX = cell(R,1); t_rawY = cell(R,1);
a_rawX = cell(R,1); a_rawY = cell(R,1);

for r=1:R
    A = Sall.(fn{r});
    if size(A,2)<2, error('%s: en az iki sütun (t, ax) olmalı.', fn{r}); end
    t0 = A(:,1); ax0 = A(:,2); ay0 = []; if size(A,2)>=3, ay0 = A(:,3); end
    [t0,iu] = unique(t0,'stable'); ax0 = ax0(iu); if ~isempty(ay0), ay0=ay0(iu); end

    [tX,ax] = regrid_to_target(t0, ax0, prep);
    if ~isempty(ay0), [~,ay] = regrid_to_target(t0, ay0, prep); else, ay=[]; end

    t_rawX{r} = tX; a_rawX{r} = ax;
    t_rawY{r} = tX; a_rawY{r} = ay;   % Y varsa aynı ızgara

    % Hedef dt kontrol uyarısı (resample_mode='off' iken dahi)
    dtX = median(diff(tX),'omitnan');
    tol = max(prep.tol_rel*max(prep.target_dt,eps), 1e-12);
    if abs(dtX - prep.target_dt) > tol
        warning('Kayıt #%d dt=%.6g s, hedef=%.6g s (resample kapalı).', r, dtX, prep.target_dt);
    end
end

%% -------------------- Şiddet eşitleme (scaled set) --------------------

t_sclX = t_rawX; t_sclY = t_rawY;
a_sclX = a_rawX; a_sclY = a_rawY;   % default: scaled = raw

if pp.on.intensity
    zeta_SA = pp.PSA.zeta; band_fac = pp.PSA.band_fac; Np_band = pp.PSA.Np_band;

    % Opsiyonel CMS hedefi
    useCMS = false; T_cms = []; Sa_cms = [];
    if pp.on.CMS && exist('cms_target.mat','file')
        Scms = load('cms_target.mat');
        if isfield(Scms,'T_cms') && isfield(Scms,'Sa_cms') && numel(Scms.T_cms)==numel(Scms.Sa_cms)
            T_cms  = Scms.T_cms(:); Sa_cms = Scms.Sa_cms(:); useCMS = true;
        end
    end

    % Önceden hesaplanan PSA var mı?
    global PSA_PREP_CACHE
    params_vec = [T1, zeta_SA, band_fac(:).', Np_band, ...
                  pp.PSA.downsample_dt, double(pp.PSA.downsample_enabled)];
    prep_key = hash_psa_prep(t_rawX, a_rawX, params_vec);

    if ~isempty(PSA_PREP_CACHE) && isfield(PSA_PREP_CACHE,'key') ...
            && strcmp(PSA_PREP_CACHE.key, prep_key)
        tPSA   = PSA_PREP_CACHE.tPSA;
        agPSA  = PSA_PREP_CACHE.agPSA;
        Sa_band = PSA_PREP_CACHE.Sa_band;
    else
        % Band Sa (parfor opsiyonel)
        tPSA   = cell(R,1);
        agPSA  = cell(R,1);
        Sa_band = zeros(R,1);
        canPar  = pp.PSA.use_parfor && ~isempty(ver('parallel'));
        if canPar
            parfor r=1:R
                if pp.PSA.downsample_enabled
                    [tLoc,agLoc] = psa_grid(t_rawX{r}, a_rawX{r}, pp.PSA.downsample_dt);
                else
                    tLoc  = t_rawX{r};
                    agLoc = a_rawX{r};
                end
                Sa_band(r) = f_band(tLoc, agLoc, T1, zeta_SA, band_fac, Np_band);
                tPSA{r}  = tLoc;
                agPSA{r} = agLoc;
            end
        else
            for r=1:R
                if pp.PSA.downsample_enabled
                    [tLoc,agLoc] = psa_grid(t_rawX{r}, a_rawX{r}, pp.PSA.downsample_dt);
                else
                    tLoc  = t_rawX{r};
                    agLoc = a_rawX{r};
                end
                Sa_band(r) = f_band(tLoc, agLoc, T1, zeta_SA, band_fac, Np_band);
                tPSA{r}  = tLoc;
                agPSA{r} = agLoc;
            end
        end
        PSA_PREP_CACHE = struct('key',prep_key,'tPSA',{tPSA}, ...
                                'agPSA',{agPSA},'Sa_band',Sa_band);
    end

    Sa_band_target = median(Sa_band);

    % Her kayıt için ölçek uygula (yalnız GA amaçlı; solver varsayılanı RAW)
    for r=1:R
        Sab_r  = f_band(tPSA{r}, agPSA{r}, T1, zeta_SA, band_fac, Np_band);
        s_band = Sa_band_target / max(Sab_r,eps);

        if useCMS
            Sa_rec = f_vec(tPSA{r}, agPSA{r}, T_cms, zeta_SA);
            nume   = max(sum(Sa_rec.*Sa_cms),eps);
            deno   = max(sum(Sa_rec.*Sa_rec),eps);
            s_cms  = nume/deno;
            s_hyb  = s_band^(1-pp.gammaCMS) * s_cms^(pp.gammaCMS);
        else
            s_hyb  = s_band;
        end

        a_sclX{r} = s_hyb * a_rawX{r};
        if ~isempty(a_rawY{r}), a_sclY{r} = s_hyb * a_rawY{r}; end
    end
end

%% -------------------- Arias pencereleri (raw & scaled) -----------------

[t5x_raw,t95x_raw,t5y_raw,t95y_raw] = deal(zeros(R,1),zeros(R,1),nan(R,1),nan(R,1));
[t5x_scl,t95x_scl,t5y_scl,t95y_scl] = deal(zeros(R,1),zeros(R,1),nan(R,1),nan(R,1));

for r=1:R
    if pp.on.arias
        [t5x_raw(r), t95x_raw(r)] = arias_win(t_rawX{r}, a_rawX{r}, 0.05, 0.95);
        if ~isempty(a_rawY{r}), [t5y_raw(r), t95y_raw(r)] = arias_win(t_rawY{r}, a_rawY{r}, 0.05, 0.95); end

        [t5x_scl(r), t95x_scl(r)] = arias_win(t_sclX{r}, a_sclX{r}, 0.05, 0.95);
        if ~isempty(a_sclY{r}), [t5y_scl(r), t95y_scl(r)] = arias_win(t_sclY{r}, a_sclY{r}, 0.05, 0.95); end
    else
        t5x_raw(r)=t_rawX{r}(1);  t95x_raw(r)=t_rawX{r}(end);
        t5x_scl(r)=t_sclX{r}(1);  t95x_scl(r)=t_sclX{r}(end);
        if ~isempty(a_rawY{r}), t5y_raw(r)=t5x_raw(r);  t95y_raw(r)=t95x_raw(r); end
        if ~isempty(a_sclY{r}), t5y_scl(r)=t5x_scl(r);  t95y_scl(r)=t95x_scl(r); end
    end
end

%% -------------------- Baz parametre setleri (parametrebulur uyumlu) ---

% Damper geometri + malzeme (BAZ)
geom.Dp    = 0.12;                 % piston çapı [m]
geom.Lgap  = 0.20;                 % etkin strok boşluğu [m]
geom.d_o   = 0.0022;               % orifis çapı [m]
geom.Lori  = 0.03;                 % orifis uzunluğu [m]
geom.Kd    = 1.8e9;                % yağın k_b (sıkışabilirlikten gelen) için ölçek [Pa]
geom.Ebody = 2.1e11;               % gövde elastisite modülü [Pa]
geom.Ap    = pi*geom.Dp^2/4;       % piston alanı [m^2] (türetilen)

% Yay (spiral) malzeme/geo (BAZ)
sh.G      = 79e9;                  % kayma modülü [Pa]
sh.d_w    = 0.018;                 % tel çapı [m]
sh.D_m    = 0.08;                  % yay ortalama çapı [m]
sh.n_turn = 26;                    % sarım sayısı [-]

% Orifis / akışkan (BAZ)  — parametrebulur ile aynı alan adları
if ~exist('orf','var')   || ~isstruct(orf),   orf   = struct(); end
if ~exist('hyd','var')   || ~isstruct(hyd),   hyd   = struct(); end
if ~exist('therm','var') || ~isstruct(therm), therm = struct(); end
if ~exist('cfg','var')   || ~isstruct(cfg),   cfg   = struct(); end
if ~exist('num','var')   || ~isstruct(num),   num   = struct(); end

% Orifis defaultları
orf.n_orf  = getfield_default(orf, 'n_orf',  2);
orf.Cd0    = getfield_default(orf, 'Cd0',    0.6);
orf.CdInf  = getfield_default(orf, 'CdInf',  0.9);
orf.Rec    = getfield_default(orf, 'Rec',    3800);
orf.p_exp  = getfield_default(orf, 'p_exp',  1.1);
orf.p_amb  = getfield_default(orf, 'p_amb',  1.0e5);
orf.cav_sf = getfield_default(orf, 'cav_sf', 1.05);

% Termal + T-bağımlı (BAZ) — parametrebulur alanları
therm.antoine_A = getfield_default(therm,'antoine_A',5.0);
therm.antoine_B = getfield_default(therm,'antoine_B',1700);
therm.antoine_C = getfield_default(therm,'antoine_C',-80);

therm.T0_C      = getfield_default(therm,'T0_C',25);
therm.Ts0_C     = getfield_default(therm,'Ts0_C',25);
therm.T_env_C   = getfield_default(therm,'T_env_C',25);

therm.hA_o_env  = getfield_default(therm,'hA_o_env',800);
therm.hA_s_env  = getfield_default(therm,'hA_s_env',600);
therm.hA_os     = getfield_default(therm,'hA_os',   800);

therm.resFactor = getfield_default(therm,'resFactor',22);  % emniyet alt sınırı
therm.cp_oil    = getfield_default(therm,'cp_oil',   1800);
therm.cp_steel  = getfield_default(therm,'cp_steel', 500);

therm.rho_ref   = getfield_default(therm,'rho_ref',  850);
therm.T_ref_C   = getfield_default(therm,'T_ref_C',   25);
therm.alpha_rho = getfield_default(therm,'alpha_rho',7e-4);
therm.beta0     = getfield_default(therm,'beta0',   1.6e9);
therm.b_beta    = getfield_default(therm,'b_beta', -0.0035);
therm.mu_ref    = getfield_default(therm,'mu_ref',   1.2);
therm.b_mu      = getfield_default(therm,'b_mu',   -0.011);

% Hidrolik / sayısı
hyd.n_parallel = getfield_default(hyd,'n_parallel', 2);
hyd.K_leak     = getfield_default(hyd,'K_leak',     0);
hyd.Lh         = getfield_default(hyd,'Lh',     0.0009);  % emniyet
hyd.Vmin_fac   = getfield_default(hyd,'Vmin_fac',0.98); % emniyet

% Numerik guardlar (parametrebulur ile aynı isimler)
num.dP_cap       = getfield_default(num,'dP_cap',       5e7);
num.mu_min_phys  = getfield_default(num,'mu_min_phys',  0.25);
num.softmin_eps  = getfield_default(num,'softmin_eps',  1e3);



%% -------------------- Türetilen parametreler --------------------------

% Alanlar
geom.Ap   = pi*geom.Dp^2/4;
orf.Ao    = orf.n_orf * (pi*geom.d_o^2/4);

% Paralel damper etkileri
nd               = max(1, getfield_default(hyd,'n_parallel',1));
geom.Ap_eff      = nd * geom.Ap;
orf.Ao_eff       = nd * orf.Ao;
hyd.n_parallel   = nd;   % garanti

% Qsat/doygunluk için referans debi (ilk kodla uyumlu)
cd_ref        = max(orf.CdInf, orf.Cd0);
dp_for_qcap   = getfield_default(num,'dP_cap', 3e8);     % num.dP_cap varsa onu kullan
Ae_ref        = max(cd_ref * orf.Ao_eff, 1e-12);
num.Qcap_big  = getfield_default(num,'Qcap_big', ...
                    hyd.Vmin_fac * Ae_ref * sqrt( 2*dp_for_qcap / max(therm.rho_ref,100) ) );


% Yay rijitlikleri
k_p   = sh.G*sh.d_w^4/(8*sh.n_turn*sh.D_m^3);   % coil spring
k_h   = geom.Kd*geom.Ap^2/geom.Lgap;            % hidrolik (tek damper)
k_s   = geom.Ebody*geom.Ap/geom.Lgap;           % gövde (tek damper)
k_hyd = 1/(1/max(k_h,eps) + 1/max(k_s,eps));    % seri birleşim
k_sd  = nd * (k_hyd + k_p);                     % kat başına efektif yay

% Yağ hacmi ve ısı kapasiteleri (parametrebulur → eq_demo_* ile uyumlu)
nStories         = n - 1;
steel_to_oil_mass_ratio = 1.5;

hyd.V0           = 0.5 * (geom.Ap * (2*geom.Lgap));        % tek hazne tahmini
V_oil_per        = therm.resFactor * (geom.Ap * (2*geom.Lgap));  % tek damper yağ [m^3]
nDtot            = nStories * nd;

m_oil_tot        = nDtot * (therm.rho_ref * V_oil_per);
m_steel_tot      = steel_to_oil_mass_ratio * m_oil_tot;
therm.C_oil      = m_oil_tot   * therm.cp_oil;    % [J/K]
therm.C_steel    = m_steel_tot * therm.cp_steel;  % [J/K]

% Bilgi: referans laminer direnç (T_ref’te)
R_lam0 = (128*therm.mu_ref*geom.Lori/(pi*geom.d_o^4)) / nd;


%% -------------------- Seçim yardımcıları ------------------------------

pickA = @(SRC,Rid,DIR) pick_series(SRC,Rid,DIR, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl);

%% ==================== GA Koşusu (isteğe bağlı) ========================

% Çok aşamalı GA: önce set-1 (geometri/rijitlik), sonra set-3 (akışkan/termal)
% Kullanım: ga.opt.multi_stage = [1 3]; // yoksa mevcut tek aşamalı akış sürer
if isfield(ga,'opt') && isfield(ga.opt,'enable') && ga.opt.enable
        % --- GA başlamadan: simulate ve yardımcılar işçilere eklensin ---
    assert(~isempty(which('simulate')), 'simulate.m path''te değil veya görünmüyor.');
pool = gcp('nocreate');
if ~isempty(pool)
    addAttachedFiles(pool, { ...
        'simulate.m', ...
        'decode_design_apply.m','ga_get_bounds.m','ensure_cfg_defaults.m', ...
        'mck_with_damper_adv.m' ...
    });
end
    % --- Aşama listesi (yoksa tek aşama: mevcut design_set) ---
    if isfield(ga.opt,'multi_stage') && ~isempty(ga.opt.multi_stage)
        stage_list = ga.opt.multi_stage(:).';
    else
        stage_list = ga.design_set;
    end

    last_best = []; last_set = NaN;

    for si = 1:numel(stage_list)
    ga.design_set = stage_list(si);


    % --- Değişken sınırları + override ---
    [lb,ub,int_idx,~] = ga_get_bounds(ga.design_set);
    if isfield(ga.opt,'lb_override') && ~isempty(ga.opt.lb_override)
        lb = max(lb, ga.opt.lb_override);
    end

    if isfield(ga.opt,'ub_override') && ~isempty(ga.opt.ub_override)
        ub = min(ub, ga.opt.ub_override);
    end
    if isempty(int_idx), IntCon = [];
    else, IntCon = int_idx(:)'; end

    % --- Bilgi satırı (banner) ---
    switch ga.design_set
      case 1
    fprintf(['[GA] set=1 | pop=%d | gen=%d | bounds(d_o)=%.1f..%.1f mm | ' ...
             'Dp=%.0f..%.0f mm | Lgap=%.0f..%.0f mm\n'], ...
        ga.opt.popsize, ga.opt.maxgen, ...
        1e3*lb(1), 1e3*ub(1), ...
        1e3*lb(8), 1e3*ub(8), ...
        1e3*lb(9), 1e3*ub(9));


      case 2
        fprintf(['[GA] set=2 | pop=%d | gen=%d | bounds(n_orf)=%g..%g | Cd0=%.2f..%.2f | CdInf=%.2f..%.2f | ' ...
                 'Rec=%g..%g | p_exp=%.1f..%.1f | cav_{sf}=%.2f..%.2f | Lh=%.1f..%.1f mm | K_{leak}=%.e..%.e | resFactor=%g..%g\n'], ...
            ga.opt.popsize, ga.opt.maxgen, ...
            lb(1),ub(1), lb(2),ub(2), lb(3),ub(3), lb(4),ub(4), lb(5),ub(5), lb(6),ub(6), ...
            1e3*lb(7),1e3*ub(7), lb(8),ub(8), lb(9),ub(9));
      case 3
        fprintf(['[GA] set=3 | pop=%d | gen=%d | bounds(n_orf)=%g..%g | Cd0=%.2f..%.2f | ' ...
                 'mu_{ref}=%.2g..%.2g Pa·s | b_{\\mu}=%.3g..%.3g 1/°C | \\beta_0=%.1e..%.1e Pa | ' ...
                 'b_{\\beta}=%.3g..%.3g 1/°C | hA_{os}=%g..%g W/K | dP_{cap}=%.1e..%.1e Pa | ' ...
                 'Vmin_{fac}=%.2f..%.2f | resFactor=%g..%g\n'], ...
            ga.opt.popsize, ga.opt.maxgen, ...
            lb(1),ub(1), lb(2),ub(2), lb(3),ub(3), lb(4),ub(4), lb(5),ub(5), lb(6),ub(6), ...
            lb(7),ub(7), lb(8),ub(8), lb(9),ub(9), lb(10),ub(10));
    end

    % --- Başlangıç nüfusu (LHS) ---
    P0 = lhs_population(lb,ub,ga.opt.popsize);
    if isfield(ga.opt,'lhs_refine') && ga.opt.lhs_refine
        P0(1,:) = 0.5*(lb+ub);
    end

    % --- Aşama-özel kısıt/PF ayarı ---
    cons_stage = cons;
    cfg_stage  = cfg;
   switch ga.design_set
    case 1  % geometri/rijitlik
        cons_stage.on.dp_quant    = true;
        cons_stage.on.thermal_dT  = true;
        cons_stage.on.cav_frac    = false;
        cons_stage.on.qsat_margin = true;

        cons_stage.pen.lambda.spring_tau     = 0.8;
        cons_stage.pen.lambda.spring_slender = 0.4;

        cfg_stage.PF.gain = 1.7;
        cfg_stage.PF.tau  = 0.9;

    case 2  % hidrolik (kavitasyonu önce bastır)
        cons_stage.on.dp_quant    = true;
        cons_stage.on.cav_frac    = true;
        cons_stage.on.qsat_margin = true;
        cons_stage.on.thermal_dT  = false;

        % cezaları biraz dengeli tut
        cons_stage.pen.lambda.cav_frac  = 1.5;
        cons_stage.pen.lambda.dp_quant  = 0.5;
        cons_stage.pen.lambda.qsat_margin = 0.3;

        % guard rails (decode sonrası yine clamp ediliyor)
        hyd.Lh          = max(hyd.Lh, 3.5e-3);
        hyd.Vmin_fac    = max(hyd.Vmin_fac, 0.93);
        therm.resFactor = max(therm.resFactor, 12);

        cfg_stage.PF.gain = 0.80;
        cfg_stage.PF.tau  = 0.9;

    case 3  % termal/akışkan tamamlayıcı
        cons_stage.on.dp_quant    = true;
        cons_stage.on.cav_frac    = true;
        cons_stage.on.qsat_margin = true;
        cons_stage.on.thermal_dT  = true;

        cons_stage.pen.lambda.thermal_dT = 0;   % aşırı ısınmayı izle ama kilitleme
        cons_stage.pen.lambda.cav_frac   = 0.2;

        % dP_cap geniş; solver güvenliği için yerinde kalsın
        num.dP_cap = max(num.dP_cap, 3e8);

        cfg_stage.PF.gain = 0.80;
        cfg_stage.PF.tau  = 0.90;
end

    % --- Fitness sarıcı ---
   inner_fitfun = @(xx) eval_fitness_for_x(xx, ga.design_set, ...
    obj, cons_stage, pp.tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg_stage, ...
    ga.opt.cache_enable, ga.opt.fail_early_k);



    fhandle = @(xx) compact_log_wrapper(xx, inner_fitfun);

    % --- GA seçenekleri ---
    opts = optimoptions('ga', ...
        'UseParallel', ga.opt.use_parallel, ...
        'InitialPopulationMatrix', P0, ...
        'PopulationSize', ga.opt.popsize, ...
        'MaxGenerations', ga.opt.maxgen, ...
        'Display','off');

    % --- Güvenli ön-tanımlar (xbest henüz yokken kullanmak YASAK) ---
    xbest = []; fbest = []; pop = []; scores = []; exitflag = [];

    % --- Çalıştır ---
    [xbest, fbest, output, pop, scores, exitflag] = ga_call_compat(fhandle, lb, ub, IntCon, opts);
    

    % --- Fallback: GA bir şey döndürmediyse orta noktayı kullan ---
    if isempty(xbest)
        xbest = 0.5*(lb+ub);
        fbest = inf;
        pop   = []; scores = [];
    end

    % --- Log & kalıcılaştır ---
    runLog = struct('stage',ga.design_set,'xbest',xbest,'fbest',fbest,'output',output, ...
                    'pop',pop,'scores',scores,'exitflag',exitflag, ...
                    'ga_options',opts,'bounds',struct('lb',lb,'ub',ub),'seed',ga.opt.seed);
    if isfield(ga.opt,'save_log') && ~isempty(ga.opt.save_log)
        try
            [pth,base,ext] = fileparts(ga.opt.save_log); if isempty(pth), pth='.'; end
            save(fullfile(pth, sprintf('%s_stage%d%s',base,ga.design_set,ext)), 'runLog');
        catch ME
            warning('runLog kaydedilemedi: %s', ME.message);
        end
    end

    % --- Bu aşamanın en iyisini uygula → bir sonraki aşama bunun üstünde arar ---
    ga.enable = true; ga.x = xbest;
    [geom, sh, orf, hyd, therm, num, ga] = decode_design_apply(ga, geom, sh, orf, hyd, therm, num);
    ga.enable = false; ga.x = [];

    last_best = xbest; last_set = ga.design_set;
    fprintf('\n=== GA Stage %d Bitti ===\nBest f = %.6g | xbest = [%s]\n', ...
        ga.design_set, fbest, join(string(xbest.'),', '));

    % --- Sınır daraltma: bir SONRAKİ aşama aynı set ise
    if si < numel(stage_list)
        next_set = stage_list(si+1);
        try
            [lb_sh, ub_sh] = shrink_bounds_from_pop(pop, scores, lb, ub, ga.opt.keep_top, ga.opt.buffer_fac);
            if ~isempty(IntCon)
                lb_sh(IntCon) = ceil(lb_sh(IntCon));
                ub_sh(IntCon) = floor(ub_sh(IntCon));
                lb_sh = max(lb, lb_sh);
                ub_sh = min(ub, ub_sh);
            end
            if next_set == ga.design_set
                ga.opt.lb_override = lb_sh;
                ga.opt.ub_override = ub_sh;
            else
                ga.opt.lb_override = [];
                ga.opt.ub_override = [];
            end
        catch
            ga.opt.lb_override = [];
            ga.opt.ub_override = [];
        end
    end

    % --- Aşama bazlı nesil sayısı (varsa) ---
    if isfield(ga.opt,'maxgen_stage') && numel(ga.opt.maxgen_stage)>=si+1
        ga.opt.maxgen = ga.opt.maxgen_stage(si+1);
    end
end


    % Son aşamanın x’i overlay/raporlar için kalsın
    ga.enable = true; ga.x = last_best; ga.best_x = last_best; ga.best_set = last_set;
    ga_dbg = struct('enable',true,'design_set',ga.best_set,'x',ga.best_x);
[geom_dbg, ~, orf_dbg, hyd_dbg, therm_dbg, ~, ~] = decode_design_apply(ga_dbg, geom, sh, orf, hyd, therm, num);

fprintf('DBG set-%d: n_orf=%d | d_o=%.3f mm | Ao=%.3e m^2 | Lgap=%.1f mm | Vmin_fac=%.2f | Lh=%.3f mm | resFactor=%.0f\n', ...
    ga.best_set, orf_dbg.n_orf, 1e3*geom_dbg.d_o, orf_dbg.n_orf*(pi*geom_dbg.d_o^2/4), ...
    1e3*geom_dbg.Lgap, hyd_dbg.Vmin_fac, 1e3*hyd_dbg.Lh, therm_dbg.resFactor);

end

% --- GA sonucunu uygula ve türetilenleri güncelle ---
if isfield(ga, 'best_x') && ~isempty(ga.best_x)
    ga.enable = true;
    ga.x = ga.best_x;
    ga.design_set = ga.best_set;
    [geom, sh, orf, hyd, therm, num, ga] = decode_design_apply(ga, geom, sh, orf, hyd, therm, num);
    ga.enable = false;
    ga.x = [];
end

% Türetilen parametreleri güncelle
geom.Ap   = pi*geom.Dp^2/4;
orf.Ao    = orf.n_orf * (pi*geom.d_o^2/4);
nd               = max(1, getfield_default(hyd,'n_parallel',1));
geom.Ap_eff      = nd * geom.Ap;
orf.Ao_eff       = nd * orf.Ao;
hyd.n_parallel   = nd;
cd_ref        = max(orf.CdInf, orf.Cd0);
dp_for_qcap   = getfield_default(num,'dP_cap', 3e8);
Ae_ref        = max(cd_ref * orf.Ao_eff, 1e-12);
num.Qcap_big  = getfield_default(num,'Qcap_big', ...
                    hyd.Vmin_fac * Ae_ref * sqrt( 2*dp_for_qcap / max(therm.rho_ref,100) ) );
k_p   = sh.G*sh.d_w^4/(8*sh.n_turn*sh.D_m^3);
k_h   = geom.Kd*geom.Ap^2/geom.Lgap;
k_s   = geom.Ebody*geom.Ap/geom.Lgap;
k_hyd = 1/(1/max(k_h,eps) + 1/max(k_s,eps));
k_sd  = nd * (k_hyd + k_p);
nStories         = n - 1;
steel_to_oil_mass_ratio = 1.5;
hyd.V0           = 0.5 * (geom.Ap * (2*geom.Lgap));
V_oil_per        = therm.resFactor * (geom.Ap * (2*geom.Lgap));
nDtot            = nStories * nd;
m_oil_tot        = nDtot * (therm.rho_ref * V_oil_per);
m_steel_tot      = steel_to_oil_mass_ratio * m_oil_tot;
therm.C_oil      = m_oil_tot   * therm.cp_oil;
therm.C_steel    = m_steel_tot * therm.cp_steel;
R_lam0 = (128*therm.mu_ref*geom.Lori/(pi*geom.d_o^4)) / nd;



%% -------------------- (1) Tek koşu: görselleştirme/diagnostic ---------

rec = min(max(1, vis.sel_rec), R);
dir = upper(string(vis.sel_dir));
[t_sim, ag_sim, t5_sim, t95_sim] = pickA(sel.sim_source, rec, dir);
% Kuyruk ekle (sim)
dt_sim  = median(diff(t_sim));
t_tail  = (t_sim(end)+dt_sim:dt_sim:t_sim(end)+pp.tail_sec).';
t_sim   = [t_sim; t_tail];
ag_sim  = [ag_sim; zeros(size(t_tail))];

% PF rampası (sim penceresine göre) — guard
cfg = set_pf_ton_if_auto(cfg, t5_sim, 0.5);

fprintf('SIM source=%s | rec #%d | dir=%s | N=%d | dt=%.3g s | Arias [%.2f, %.2f] s\n', ...
        sel.sim_source, rec, dir, numel(t_sim), dt_sim, t5_sim, t95_sim);

[x0_sim,a0_sim] = lin_MCK_consistent(t_sim, ag_sim, M, Cstr, K);
[xD_sim, aD_sim, dlog_sim, vD_sim] = mck_with_damper_adv( ...
    t_sim, ag_sim, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg);
% Grafik seti (opsiyonel ikinci koşu)
use_plot_run = false;
cfg_plot = cfg;   % tern(...) çağrısı için hazır dursun
if vis.dual_run_if_needed && ~strcmpi(sel.plot_source, sel.sim_source)
    [t_plot, ag_plot, t5_plot, t95_plot] = pickA(sel.plot_source, rec, dir);
    dt_plot = median(diff(t_plot));
    t_tail  = (t_plot(end)+dt_plot:dt_plot:t_plot(end)+pp.tail_sec).';
    t_plot  = [t_plot; t_tail];
    ag_plot = [ag_plot; zeros(size(t_tail))];

    cfg_plot = set_pf_ton_if_auto(cfg_plot, t5_plot, 0.5);

    [x0,a0] = lin_MCK_consistent(t_plot, ag_plot, M, Cstr, K);
    [xD,aD,dlog,vD] = mck_with_damper_adv( ...
        t_plot, ag_plot, M, Cstr, K, k_sd, geom, orf, hyd, therm, num, cfg_plot);

    use_plot_run = true;
    fprintf('PLOT source=%s | rec #%d | dir=%s | N=%d | dt=%.3g s | Arias [%.2f, %.2f] s\n', ...
            sel.plot_source, rec, dir, numel(t_plot), dt_plot, t5_plot, t95_plot);
else
    % Plot, sim ile aynı
t_plot=t_sim; ag_plot=ag_sim; t5_plot=t5_sim; t95_plot=t95_sim;
x0=x0_sim; a0=a0_sim; xD=xD_sim; aD=aD_sim; dlog=dlog_sim; vD=vD_sim;
dt_plot = dt_sim;   % başlıktaki dt için

end

% ---- Enerji/diagnostik (sadece yazdırma) ----
active_cfg  = tern(use_plot_run, cfg_plot, cfg);
dlog_active = tern(use_plot_run, dlog, dlog_sim);

E_orif  = dlog_active.E_cum(end);                 % orifis enerjisi (J)
P_struc = sum( (vD * Cstr) .* vD, 2 );
E_struc = trapz(t_plot, P_struc);

fprintf('E_orifice = %.3e J | E_struct = %.3e J | oran = %.3f\n', ...
        E_orif, E_struc, E_orif/max(E_struc,eps));

P_mech = sum( (dlog_active.F_story .* (vD(:,2:end)-vD(:,1:end-1))), 2 );
fprintf('⟨P_mech⟩ = %.3e W (negatif ise net sönüm)\n', mean(P_mech,'omitnan'));

fprintf('CHECK → leak=%.2e, Lh=%.3e, Ao=%.2e m^2, PF=%s\n', ...
    hyd.K_leak, hyd.Lh, orf.n_orf*(pi*geom.d_o^2/4), active_cfg.PF.mode);

fprintf('dP_orf q95 = %.3e Pa | Q95 = %.3e m^3/s\n', ...
    prctile(dlog_active.dP_orf_time_max,95), dlog_active.Q_abs_p95);

if isfield(dlog_active,'cav_margin_min')
    fprintf('cav_margin_min = %.1f kPa\n', 1e-3*dlog_active.cav_margin_min);
end

% === Tek FIGÜR (ilk koddaki stil): üstte yer değiştirme, altta mutlak ivme ===
if vis.make_plots
    j_mon = min(max(1, obj.idx_disp_story), size(xD,2));   % izlenen kat

    % Mutlak ivme (dampersiz/damperli)
    a_abs0 = a0 + ag_plot * ones(1, size(a0,2));
    a_absD = aD + ag_plot * ones(1, size(aD,2));

    figure('Color','w','Visible','on','Name', ...
        sprintf('%d. kat: dampersiz vs damperli (TERMAL+Cd(Re)+Kaçak+Vmin+NUM+β(T)+P_{sat})', j_mon));

    % --- Üst: yer değiştirme ---
    subplot(2,1,1); hold on; grid on;
    plot(t_plot, x0(:, j_mon), 'k--', 'DisplayName','dampersiz');
    plot(t_plot, xD(:, j_mon), 'b-',  'DisplayName','damperli');
    ylabel(sprintf('x_{%d} [m]', j_mon)); legend('show','Location','best');
    title(sprintf('N=%d, dt=%.4fs | d_o=%.1f mm, gain=%.2f, \\tau=%.2f s', ...
        numel(t_plot), dt_plot, 1e3*geom.d_o, active_cfg.PF.gain, active_cfg.PF.tau));

    % --- Alt: mutlak ivme ---
    subplot(2,1,2); hold on; grid on;
    plot(t_plot, a_abs0(:, j_mon), 'k--', 'DisplayName','dampersiz');
    plot(t_plot, a_absD(:, j_mon), 'b-',  'DisplayName','damperli');
    ylabel(sprintf('a_{%d} [m/s^2]', j_mon)); xlabel('t [s]'); legend('show','Location','best');
end

%% -------------------- (2) Amaç fonksiyonu — kompakt (ilk koddaki gibi) ----

% Girdi olarak şunlar zaten mevcut olmalı:
% t_plot, dt_plot, vD, dlog_active, M, K, k_sd, obj, cons, orf, hyd, therm, num, t5_plot, t95_plot

% Kapasiteler
dPcap_eff = (isfield(num,'dP_cap') && isfinite(num.dP_cap)) * num.dP_cap + ...
            (~(isfield(num,'dP_cap') && isfinite(num.dP_cap))) * 3e8;

cd_ref   = max(orf.CdInf, orf.Cd0);
Ae_ref   = cd_ref * max(orf.Ao_eff, 1e-12);                 % paralel ve Cd dahil
rho_ref  = max(therm.rho_ref, 100);
Qcap_ref = hyd.Vmin_fac * Ae_ref * sqrt(2*dPcap_eff / rho_ref);

% Gözlenenler (dlog_active üzerinden)
Qp95_max = dlog_active.Q_abs_p95;                            % m^3/s
dp95_max = prctile(dlog_active.dP_orf_time_max, 95);         % Pa

% Uygunluk bayrakları
Qmargin = (isfield(cons,'hyd') && isfield(cons.hyd,'Q_margin') && isfinite(cons.hyd.Q_margin)) ...
            * cons.hyd.Q_margin + (~(isfield(cons,'hyd') && isfield(cons.hyd,'Q_margin'))) * 0.90;
okQ  = (Qp95_max <= Qmargin * Qcap_ref);
okdp = (dp95_max <= dPcap_eff);

% zeta_eq tahmini (1. mod)
[PHI,LAM] = eig(K,M);
WW = sqrt(builtin('diag', LAM));               % 'diag' fonksiyonu gölgelenmesin diye builtin

w1   = min(WW);
phi1 = PHI(:, WW==w1); phi1 = phi1(:,1);
phi1 = phi1 / max(abs(phi1));
m_eff = phi1.'*M*phi1;
k_eff = phi1.'*K*phi1;

% Arias penceresinde izlenen katın hız RMS'i
j_mon = min(max(1, obj.idx_disp_story), size(vD,2));
maskA = (t_plot >= t5_plot) & (t_plot <= t95_plot);
v_rms = sqrt(mean(vD(maskA, j_mon).^2, 'omitnan'));
Xeq   = v_rms / max(w1, 1e-6);
Est   = 0.5 * (k_eff + k_sd) * Xeq^2;

% Mekanik güçten efektif sönüm (yalnız sönümleyici kısım)
P_mech  = sum( dlog_active.F_story .* (vD(:,2:end) - vD(:,1:end-1)), 2 );
P_diss  = max(-P_mech, 0);                                  % negatif → sönüm kabulü
zeta_eq = sum(P_diss) * dt_plot / max(4*pi*Est, eps);

% <P_mech> işareti (pozitif ⇒ net sönüm)
Pmech_avg = -mean(P_mech, 'omitnan');
okP = (Pmech_avg > 0);

% Skor
score = zeta_eq ...
        - 0.7*double(~okP) ...
        - 0.5*double(~okQ) ...
        - 0.5*double(~okdp);

fprintf('\n=== AMAÇ (kompakt) ===\n');
fprintf('zeta_eq=%.3f | <P_mech>_diss=%.3e W | Q95=%.3e | dp95=%.3e | ok=[Q %d, dp %d, P %d] | score=%.3f\n', ...
    zeta_eq, Pmech_avg, Qp95_max, dp95_max, okQ, okdp, okP, score);

% GA minimizasyonu ile uyum için amaç değeri
J = -score;

%% -------------------- (3) Kısıtlar: Penalty hesap (R kayıt) ----------
[Penalty, cons_detail] = evaluate_constraints_over_records( ...
    cons, cons.src_for_constraints, obj, pp.tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, tern(ga.enable,ga.design_set,0), tern(ga.enable,ga.x,[]));

Fitness = J + Penalty;
fprintf('λpen = [tau=%.2f, stroke=%.2f, dT=%.2f, cav=%.2f]\n', ...
    cons.pen.lambda.spring_tau, cons.pen.lambda.stroke, ...
    cons.pen.lambda.thermal_dT, cons.pen.lambda.cav_frac);
fprintf('PF: mode=%s, t_on=%.2fs, tau=%.2f, gain=%.2f\n', ...
    cfg.PF.mode, cfg.PF.t_on, cfg.PF.tau, cfg.PF.gain);



fprintf('\n================ Kısıt & Ceza Sonuçları ================\n');
fprintf('Penalty = %.6g  |  Fitness = J + Penalty = %.6g\n', Penalty, Fitness);
% === Kısıt/cihaz özeti → CSV ===
try, mkdir('out'); end

ratio = cons_detail.ratios;
okflag = @(r) (r<=1+1e-12);

Names   = {'spring_tau','spring_slender','stroke','force_cap','dp_quant','thermal_dT','cav_frac','qsat_margin'};
Limits  = [cons.spring.tau_allow, cons.spring.lambda_max, cons.stroke.util_factor*geom.Lgap, ...
           cons.force.F_cap, num.dP_cap, cons.thermal.cap_C, cons.hyd.cav_frac_cap, cons.hyd.Q_margin * getfield_default(num,'Qcap_big', ...
           0.4 * ( max(max(orf.CdInf,orf.Cd0)*orf.Ao_eff, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) )];

Ratios  = [ratio.spring_tau, ratio.spring_slender, ratio.stroke, ratio.force_cap, ...
           ratio.dp_quant, ratio.thermal_dT, ratio.cav_frac, ratio.qsat_margin];

Flags = strings(numel(Ratios),1);
okmask = Ratios <= (1 + 1e-12);
Flags(okmask)  = "OK";
Flags(~okmask) = "VIOL";


Vals_Fmax   = max(cons_detail.Fmax_records,   [], 'omitnan');
Vals_stroke = max(cons_detail.stroke_records, [], 'omitnan');
Vals_dpq    = max(cons_detail.dpq_records,    [], 'omitnan');
Vals_dT     = max(cons_detail.dT_records,     [], 'omitnan');
Vals_cav    = max(cons_detail.cav_records,    [], 'omitnan');
Vals_Qp95   = max(cons_detail.Qp95_records,   [], 'omitnan');

T = table(Names.', Ratios.', Flags, ...
    'VariableNames', {'Constraint','Ratio','Flag'});

% Ek cihaz metrikleri ayrı tablo: (zarf değerleri)
T_dev = table(Vals_Fmax, Vals_stroke, Vals_dpq, Vals_dT, Vals_cav, Vals_Qp95, ...
    'VariableNames', {'Fmax_N','stroke_m','dp_q95_Pa','dT_C','cav95','Q95_m3s'});

writetable(T,    fullfile('out','cons_summary.csv'));
writetable(T_dev,fullfile('out','device_summary.csv'));
fprintf('CSV yazıldı: out/cons_summary.csv, out/device_summary.csv\n');


% ---------- Yardımcılar ----------
hinge  = @(r) max(0, r - 1);
norm0  = @(x) (abs(x) < 1e-12) * 0 + (abs(x) >= 1e-12) .* x;  % -0.000 yerine 0.000
pwr    = cons.pen.power;
lam    = cons.pen.lambda;

pen_sum = 0;   % bileşenlerden yeniden hesaplanan toplam (kontrol amaçlı)

% ---------- Özet oranlar + bireysel ceza katkıları ----------
cdt = cons_detail;  % kısaltma

if cons.on.spring_tau
    r = norm0(cdt.ratios.spring_tau);
    pen_i = lam.spring_tau * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('τ_max/τ_allow         = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.spring_slender
    r = norm0(cdt.ratios.spring_slender);
    pen_i = lam.spring_slender * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('L_free/D_m / λ_max    = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.stroke
    r = norm0(cdt.ratios.stroke);
    pen_i = lam.stroke * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('stroke/(0.9*L_gap)    = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.force_cap
    r = norm0(cdt.ratios.force_cap);
    pen_i = lam.force_cap * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('F_max/F_cap           = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.dp_quant
    r = norm0(cdt.ratios.dp_quant);
    pen_i = lam.dp_quant * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('q_Δp/dP_cap           = %.3f (q=%.3f)   [pen=%.3g]\n', r, cons.dp.q, pen_i);
end

if cons.on.thermal_dT
    r = norm0(cdt.ratios.thermal_dT);
    pen_i = lam.thermal_dT * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('ΔT/ΔT_cap             = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.cav_frac
    r = norm0(cdt.ratios.cav_frac);
    pen_i = lam.cav_frac * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('cav95/γ_cap           = %.3f   [pen=%.3g]\n', r, pen_i);
end

if cons.on.qsat_margin
    r = norm0(cdt.ratios.qsat_margin);
    pen_i = lam.qsat_margin * hinge(r)^pwr; pen_sum = pen_sum + pen_i;
    fprintf('Q95/(margin*Qcap_big) = %.3f   [pen=%.3g]\n', r, pen_i);
end


% BigM (sim başarısızlığı) bilgilendirmesi ve katkısı
if cons.on.fail_bigM && cdt.any_fail
    pen_BigM = cons.pen.bigM * lam.fail_bigM;
    pen_sum  = pen_sum + pen_BigM;
    fprintf('Sim başarısız (≥1 koşu): BigM cezası eklendi. [pen=%.3g]\n', pen_BigM);
end

fprintf('--- Penalty (bileşenlerden) ≈ %.6g\n', pen_sum);



% ---------- İsteğe bağlı tanılama: zarf/metrik özetleri ----------
if isfield(cons_detail,'Fmax_records')
    % Kayıt zarfı / agregeler (dpq_all için agresyon kuralını uygula)
    Fmax_all   = max(cons_detail.Fmax_records,   [], 'omitnan');
    stroke_all = max(cons_detail.stroke_records, [], 'omitnan');
    if isfield(cons,'dp') && isfield(cons.dp,'agg') && strcmpi(cons.dp.agg,'cvar')
        dpq_all = cvar_from_samples(cons_detail.dpq_records(:), cons.alpha_CVaR_cons);
    else
        dpq_all = max(cons_detail.dpq_records, [], 'omitnan');
    end
    dT_all   = max(cons_detail.dT_records,   [], 'omitnan');
    cav_all  = max(cons_detail.cav_records,  [], 'omitnan');
    Qp95_all = max(cons_detail.Qp95_records, [], 'omitnan');

    fprintf('--- Tanılama (zarf değerleri) ---\n');
    fprintf('Fmax=%.3e N | stroke=%.3e m | dpq=%.3e Pa | ΔT=%.2f C | cav95=%.3f | Q95=%.3e m^3/s\n', ...
        Fmax_all, stroke_all, dpq_all, dT_all, cav_all, Qp95_all);
end


% ======================================================================
%                             FONKSİYONLAR
% ======================================================================
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

function y = tern(cond,a,b)
    if cond
        y = a;
    else
        y = b;
    end
end


function [J, out] = compute_objective_over_records( ...
src, obj, tail_sec, ...
t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    if numel(varargin) >= 2
        design_set = varargin{1}; x_ga = varargin{2};
    else
        design_set = 0; x_ga = [];
    end
    ...
    % local_design_ratios_one_dir içinde:
    % Simülasyon penceresi: aynı veriye kuyruk ekleyelim (PF ramp guard uyumu)
    % Kaynak seçici
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R = numel(tX);
    out_rec(R) = struct('d_rel',NaN,'a_rel',NaN,'J_r',NaN);
    logs = {};
    fail_ratio = 5;  % fallback multiplier when simulation fails

    % Önce REFERANSLAR: (damper yok, μ=1.0) → her kayıt & (X,Y mevcutsa) yön için sabitlenir
    ref = struct('X',struct('d',nan(R,1),'a',nan(R,1)), ...
                 'Y',struct('d',nan(R,1),'a',nan(R,1)));

    for r=1:R
        % X yönü referans
        if ~isempty(aX{r})
            [dref, aref] = local_ref_metrics(tX{r}, aX{r}, t5x(r), t95x(r), obj, M,Cstr,K, n);
            ref.X.d(r)=dref; ref.X.a(r)=aref;
        end
        % Y yönü referans (varsa)
        if ~isempty(aY{r})
            [dref, aref] = local_ref_metrics(tY{r}, aY{r}, t5y(r), t95y(r), obj, M,Cstr,K, n);
            ref.Y.d(r)=dref; ref.Y.a(r)=aref;
        end
    end

    % Sonra TASARIM: her kayıt → yön zarfı → μ-agg → J_r
    for r=1:R
        [d_rel_X, a_rel_X] = local_design_ratios_one_dir('X', r);
        [d_rel_Y, a_rel_Y] = local_design_ratios_one_dir('Y', r);

        switch lower(obj.dir_mode)
            case 'xonly', d_rel = d_rel_X; a_rel = a_rel_X;
            case 'yonly', d_rel = d_rel_Y; a_rel = a_rel_Y;
            otherwise     % 'envelope'
                d_rel = max([d_rel_X, d_rel_Y], [], 'omitnan');
                a_rel = max([a_rel_X, a_rel_Y], [], 'omitnan');
        end

        % Ağırlıklı toplam
        J_r = obj.weights_da(1)*d_rel + obj.weights_da(2)*a_rel;

        out_rec(r).d_rel = d_rel;
        out_rec(r).a_rel = a_rel;
        out_rec(r).J_r   = J_r;
    end

    % CVaR(α) hesap
    alpha = min(max(obj.alpha_CVaR,eps),0.99);
    Jlist = [out_rec.J_r].';
    J = cvar_from_samples(Jlist, alpha);
    if ~isfinite(J)
        msg = 'Objective J is NaN/Inf';
        logs{end+1} = msg; log_msg('error', msg);
        J = 1e9;
    end

    out = struct('records', out_rec, 'logs', {logs});

    % ---- iç yardımcılar ----
 
    function [d_ref, a_ref] = local_ref_metrics(t, ag, t5, t95, obj, M,C,K, n)
        % Kuyruk eklemeden baz referans (damper yok)
        [x0,a0] = lin_MCK_consistent(t, ag, M, C, K);
        d_ref = max(abs(x0(:,min(obj.idx_disp_story,n))));
        a_ref = acc_metric_from_series(t, a0(:,min(obj.idx_acc_story,n)), t5, t95, obj);
        d_ref = max(d_ref, eps);
        a_ref = max(a_ref, eps);
    end

    function [d_rel_dir, a_rel_dir] = local_design_ratios_one_dir(which, r)
       
    switch upper(string(which))
        case "X"
            t=tX{r}; ag=aX{r}; t5=t5x(r); t95=t95x(r);
            d_ref=ref.X.d(r); a_ref=ref.X.a(r);
        case "Y"
            t=tY{r}; ag=aY{r}; t5=t5y(r); t95=t95y(r);
            d_ref=ref.Y.d(r); a_ref=ref.Y.a(r);
    end
    if isempty(ag) || ~isfinite(d_ref) || ~isfinite(a_ref)
        d_rel_dir = NaN; a_rel_dir = NaN; return;
    end


 
    % --- Simülasyon penceresi: aynı veriye kuyruk ekleyelim ---
 
    dt = median(diff(t));
tail_sec_loc = tail_sec;
    t_tail = (t(end)+dt : dt : t(end)+tail_sec_loc).';
    t_s    = [t;  t_tail];
    ag_s   = [ag; zeros(size(t_tail))];

    % PF ramp t_on: Arias t5 + 3
    cfg_dir = set_pf_ton_if_auto(cfg, t5, 0.5);


    % μ senaryoları
    mus   = obj.mu_scenarios(:).';
    d_vals = nan(size(mus));  a_vals = nan(size(mus));

    for k = 1:numel(mus)
        resp = simulate( ...
            design_set, x_ga, mus(k), t_s, ag_s, ...
            M, Cstr, K, n, geom, sh, orf, hyd, therm, num, cfg_dir);

        if ~resp.ok
    msg = sprintf('simulate fail (rec=%d dir=%s mu=%g): %s', r, which, mus(k), resp.msg);
    logs{end+1} = msg; %#ok<AGROW>
    log_msg('warn', msg);
    d_vals(k) = fail_ratio * d_ref;
    a_vals(k) = fail_ratio * a_ref;
    continue;
end

        x  = resp.y.x(1:numel(t), :);   % kuyruğu at
        aR = resp.y.a(1:numel(t), :);

        d  = max(abs(x(:,min(obj.idx_disp_story,n))));
        am = acc_metric_from_series(t, aR(:,min(obj.idx_acc_story,n)), t5, t95, obj);

        d_vals(k) = d;  a_vals(k) = am;
    end

    % μ-aggregation
    switch lower(obj.mu_aggregate)
        case 'weighted'
            w = obj.mu_weights(:); w = w/sum(w);
            d_agg = nansum(w.*d_vals(:));
            a_agg = nansum(w.*a_vals(:));
        otherwise
            d_agg = max(d_vals, [], 'omitnan');
            a_agg = max(a_vals, [], 'omitnan');
    end

    d_rel_dir = d_agg / max(d_ref, eps);
    a_rel_dir = a_agg / max(a_ref, eps);
    if ~isfinite(d_rel_dir)
        d_rel_dir = fail_ratio;
        msg = sprintf('NaN d_rel at rec=%d dir=%s', r, which);
        logs{end+1} = msg; log_msg('warn', msg);
    end
    if ~isfinite(a_rel_dir)
        a_rel_dir = fail_ratio;
        msg = sprintf('NaN a_rel at rec=%d dir=%s', r, which);
        logs{end+1} = msg; log_msg('warn', msg);
    end
end

end

function val = acc_metric_from_series(t, a, t5, t95, obj)
    % Arias penceresi içi metrikler
    if obj.use_arias_window
        w = (t>=t5 & t<=t95);
    else
        w = true(size(t));
    end
    aa = a(w); tt = t(w);
    if isempty(aa) || numel(aa)<2
        val = NaN; return;
    end
    switch lower(obj.acc_metric)
        case 'energy'
            val = trapz(tt, aa.^2);                 % enerji
        case 'rms+p95'
            rmsv = sqrt(mean(aa.^2,'omitnan'));
            p95  = prctile(abs(aa),95);
            val  = rmsv + obj.p95_penalty_w * p95;  % hibrit
        otherwise % 'rms'
            val = sqrt(mean(aa.^2,'omitnan'));
    end
end

function v = cvar_from_samples(x, alpha)
    % Örneklerden CVaR_α (Average Value-at-Risk)
    x = x(:); x = x(isfinite(x));
    if isempty(x), v = NaN; return; end
    q = quantile(x, 1-alpha);
    tail = x(x>=q);  % büyük-kötü kuyruk
    if isempty(tail), v = q; else, v = mean(tail); end
end
function [J1, out] = compute_J1_IDR_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    % Optional GA args
    if numel(varargin) >= 2, design_set = varargin{1}; x_ga = varargin{2};
    else, design_set = 0; x_ga = []; end

    % Kaynak seçimi
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R  = numel(tX);
    Ns = n - 1;
    if isscalar(h_story_m), h_story = repmat(h_story_m, Ns, 1);
    else,                   h_story = h_story_m(:); end

    out = struct('d_rel', nan(R,1));
    mus = obj.mu_scenarios(:).';
    isWeighted = strcmpi(obj.mu_aggregate,'weighted');
    if isWeighted
        wmu = obj.mu_weights(:); wmu = wmu / max(sum(wmu), eps);
    end

    for r = 1:R
        % --- X yönü referans IDR
        dref_X = NaN; dagg_X = NaN;
        if ~isempty(aX{r})
            [dref_X, dagg_X] = local_idr_ref_and_agg(tX{r}, aX{r}, t5x(r), t95x(r), 'X', r);
        end

        % --- Y yönü (varsa)
        dref_Y = NaN; dagg_Y = NaN;
        if ~isempty(aY{r})
            [dref_Y, dagg_Y] = local_idr_ref_and_agg(tY{r}, aY{r}, t5y(r), t95y(r), 'Y', r);
        end

        % --- Yön zarfı
        switch lower(obj.dir_mode)
            case 'xonly', d_ref = dref_X; d_agg = dagg_X;
            case 'yonly', d_ref = dref_Y; d_agg = dagg_Y;
            otherwise
                d_ref = max([dref_X, dref_Y], [], 'omitnan');
                d_agg = max([dagg_X, dagg_Y], [], 'omitnan');
        end

        out.d_rel(r) = d_agg / max(d_ref, eps);
    end

    J1 = cvar_from_samples(out.d_rel, obj.alpha_CVaR);

    % ---- local helper
    function [d_ref, d_agg] = local_idr_ref_and_agg(t, ag, t5, t95, dirStr, rid)
        % Referans (damper yok, μ=1.0), pencere içinde IDR zarfı
        [x0,~] = lin_MCK_consistent(t, ag, M, Cstr, K);
        w = (t >= t5 & t <= t95);
        drift0 = x0(:,2:end) - x0(:,1:end-1);                       % Nt x Ns
        idr0   = abs(drift0(w,:)) ./ (ones(sum(w),1) * h_story(:)'); % Nt_w x Ns
        d_ref  = max(idr0,[],'all');                                 % skaler

        % Tasarım: kuyruk ekle + PF ramp guard
        dt     = median(diff(t));
        t_tail = (t(end)+dt : dt : t(end)+tail_sec).';
        t_s    = [t; t_tail];
        ag_s   = [ag; zeros(size(t_tail))];

        cfg_dir = set_pf_ton_if_auto(cfg, t5, 0.5);

        vals = nan(size(mus));
        for k = 1:numel(mus)
            resp = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                            M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir);
            if ~resp.ok
                vals(k) = 5 * d_ref;  % fail durumunda güvenli büyük ceza
                continue;
            end
            xD   = resp.y.x(1:numel(t), :);     % kuyruğu at
            driftD = xD(:,2:end) - xD(:,1:end-1);
            idrD   = abs(driftD(w,:)) ./ (ones(sum(w),1) * h_story(:)');
            vals(k)= max(idrD,[],'all');
        end

        if isWeighted, d_agg = nansum(wmu .* vals(:));
        else,          d_agg = max(vals, [], 'omitnan'); end
    end
end


function [J2, out] = compute_J2_ACC_over_records( ...
    src, obj, h_story_m, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    if numel(varargin) >= 2, design_set = varargin{1}; x_ga = varargin{2};
    else, design_set = 0; x_ga = []; end
     
 
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R   = numel(tX);
    w_p = obj.p95_penalty_w;
    out = struct('a_rel', nan(R,1));

    mus = obj.mu_scenarios(:).';
    isWeighted = strcmpi(obj.mu_aggregate,'weighted');
    if isWeighted
        wmu = obj.mu_weights(:); wmu = wmu / max(sum(wmu), eps);
    end

    for r = 1:R
        % --- X yönü: A_env ref ve agg
        aref_X = NaN; aagg_X = NaN;
        if ~isempty(aX{r})
            [aref_X, aagg_X] = local_acc_ref_and_agg(tX{r}, aX{r}, t5x(r), t95x(r), 'X', r);
        end

        % --- Y yönü (varsa)
        aref_Y = NaN; aagg_Y = NaN;
        if ~isempty(aY{r})
            [aref_Y, aagg_Y] = local_acc_ref_and_agg(tY{r}, aY{r}, t5y(r), t95y(r), 'Y', r);
        end

        % --- Yön zarfı
        switch lower(obj.dir_mode)
            case 'xonly', a_ref = aref_X; a_agg = aagg_X;
            case 'yonly', a_ref = aref_Y; a_agg = aagg_Y;
            otherwise
                a_ref = max([aref_X, aref_Y], [], 'omitnan');
                a_agg = max([aagg_X, aagg_Y], [], 'omitnan');
        end

        out.a_rel(r) = a_agg / max(a_ref, eps);
    end

    J2 = cvar_from_samples(out.a_rel, obj.alpha_CVaR);

    % ---- local helper
    function [A_ref, A_agg] = local_acc_ref_and_agg(t, ag, t5, t95, dirStr, rid)
        % Referans (damper yok) mutlak ivme zarfı (RMS+p95)
        [~,a_rel0] = lin_MCK_consistent(t, ag, M, Cstr, K);  % relatif
        a_abs0 = a_rel0 + ag(:) * ones(1, size(a_rel0,2));  % her kat için mutlak
        w = (t >= t5 & t <= t95);
        A_i = zeros(1,size(a_abs0,2));
        for i=1:size(a_abs0,2)
            ai = a_abs0(w,i);
            rmsv = sqrt(mean(ai.^2,'omitnan'));
            p95  = prctile(abs(ai),95);
            A_i(i) = rmsv + w_p * p95;
        end
        A_ref = max(A_i);   % kat zarfı

        % Tasarım: kuyruk ekle + PF guard
        dt     = median(diff(t));
        t_tail = (t(end)+dt : dt : t(end)+tail_sec).';
        t_s    = [t; t_tail];
        ag_s   = [ag; zeros(size(t_tail))];
        cfg_dir = set_pf_ton_if_auto(cfg, t5, 0.5);

        vals = nan(size(mus));
        for k = 1:numel(mus)
            resp = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                            M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir);
            if ~resp.ok
                vals(k) = 5 * A_ref;  % fail durumunda güvenli büyük ceza
                continue;
            end
            a_relD = resp.y.a(1:numel(t), :);           % relatif, kuyruğu at
            a_absD = a_relD + ag(:) * ones(1,size(a_relD,2));

            Ai = zeros(1,size(a_absD,2));
            for i=1:size(a_absD,2)
                ai = a_absD(w,i);
                rmsv = sqrt(mean(ai.^2,'omitnan'));
                p95  = prctile(abs(ai),95);
                Ai(i)= rmsv + w_p * p95;
            end
            vals(k) = max(Ai);   % kat zarfı
        end

        if isWeighted, A_agg = nansum(wmu .* vals(:));
        else,          A_agg = max(vals, [], 'omitnan'); end
    end
end

function [Penalty, out] = evaluate_constraints_over_records( ...
    cons, src, obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, varargin)

    if numel(varargin) >= 2
        design_set = varargin{1}; x_ga = varargin{2};
    else
        design_set = 0; x_ga = [];
    end
    ...

      % Kaynak seçimi
    useScaled = strcmpi(src,'scaled');
    if useScaled
        tX=t_sclX; tY=t_sclY; aX=a_sclX; aY=a_sclY;
        t5x=t5x_scl; t95x=t95x_scl; t5y=t5y_scl; t95y=t95y_scl;
    else
        tX=t_rawX;  tY=t_rawY;  aX=a_rawX;  aY=a_rawY;
        t5x=t5x_raw; t95x=t95x_raw; t5y=t5y_raw; t95y=t95y_raw;
    end

    R   = numel(tX);
    mus = cons.mu_scenarios(:).';
    logs = {};

    % ----- Erken çıkış sayaç ayarları -----
    fail_early_k = 0;
    if isfield(cons,'fail_early_k') && ~isempty(cons.fail_early_k) && isfinite(cons.fail_early_k)
        fail_early_k = max(0, round(cons.fail_early_k));
    end
    fail_count_global = 0;
    early_break_all   = false;

    % Toplayıcılar
    Fmax_records   = nan(R,1);
    stroke_records = nan(R,1);
    dpq_records    = nan(R,1);
    dT_records     = nan(R,1);
    cav_records    = nan(R,1);
    Qp95_records   = nan(R,1);
    any_fail_rec   = false(R,1);

    for r = 1:R
        % ---- Yön zarfı için X/Y ölçüleri
        rec_metrics = struct('Fmax',[],'stroke',[],'dpq',[],'dT',[],'cav95',[],'Qp95',[],'fail',[]);

        for dir = ["X","Y"]
            % Y yönü yoksa atla
            if dir=="Y" && (isempty(aY{r}) || all(isnan(aY{r})))
                continue;
            end

            % Kayıt/yön serileri
            if dir=="X"
                t = tX{r}; ag = aX{r}; t5=t5x(r); t95=t95x(r);
            else
                t = tY{r}; ag = aY{r}; t5=t5y(r); t95=t95y(r);
            end

            % ---- Tail ekle (simülasyon penceresi)
            dt    = median(diff(t));
            t_tail= (t(end)+dt : dt : t(end)+tail_sec).';
            t_s   = [t; t_tail];
            ag_s  = [ag; zeros(size(t_tail))];

            % ---- PF ramp koruması
cfg_dir = set_pf_ton_if_auto(cfg, t5, 0.5);


            % ---- μ-senaryosu zarfı
            Fmax_mu   = -Inf; stroke_mu = -Inf; dpq_mu = -Inf; dT_mu = -Inf; cav_mu = -Inf; Qp95_mu = -Inf;
            fail_mu   = false;

            % >>>>>>>>>>>>>>>>> μ DÖNGÜSÜ (BURADA) <<<<<<<<<<<<<<<<<
            for k = 1:numel(mus)
                resp = simulate(design_set, x_ga, mus(k), t_s, ag_s, ...
                                M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg_dir);

                if ~resp.ok
                    fail_mu = true;
                    fail_count_global = fail_count_global + 1;
                    msg = sprintf('simulate fail (rec=%d dir=%s mu=%g): %s', r, dir, mus(k), resp.msg);
                    logs{end+1} = msg; %#ok<AGROW>
                    log_msg('warn', msg);

                    if fail_early_k > 0 && fail_count_global >= fail_early_k
                        early_break_all = true;   % tüm döngülerden çık
                        break;                    % μ-döngüsünü kır
                    end
                    continue;  % sıradaki μ
                end

                % ---- ölçümler (başarılı koşu)
                Fmax_mu   = max(Fmax_mu,   resp.F_max);
                stroke_mu = max(stroke_mu, resp.stroke_max);
                dpq_mu    = max(dpq_mu,    resp.dP_q_time(cons.dp.q));
                dT_mu     = max(dT_mu,     resp.dT_est);
                mask = (t_s >= t5) & (t_s <= t95);
    if any(mask)
        cav_p95_win = prctile(resp.diag.cav_frac_t(mask), 95);
    else
        cav_p95_win = prctile(resp.diag.cav_frac_t, 95);  % emniyetli geri dönüş
    end
    cav_mu = max(cav_mu, cav_p95_win);
                Qp95_mu   = max(Qp95_mu,   resp.metrics.Q_abs_p95);
            end
            
            % >>>>>>>>>>>>>>>>> μ DÖNGÜSÜ BİTTİ <<<<<<<<<<<<<<<<<<<<

            % μ-döngüsü eşik nedeniyle kırıldıysa, yön döngüsünü de bırak
            if early_break_all
                break;
            end

            % Bu yönün metriklerini biriktir
            rec_metrics.Fmax   = [rec_metrics.Fmax,   Fmax_mu];
            rec_metrics.stroke = [rec_metrics.stroke, stroke_mu];
            rec_metrics.dpq    = [rec_metrics.dpq,    dpq_mu];
            rec_metrics.dT     = [rec_metrics.dT,     dT_mu];
            rec_metrics.cav95  = [rec_metrics.cav95,  cav_mu];
            rec_metrics.Qp95   = [rec_metrics.Qp95,   Qp95_mu];
            rec_metrics.fail   = [rec_metrics.fail,   fail_mu];
        end

        % Yön döngüsünden erken çıktıysak kalan kayıtları INF kabul et
        if early_break_all
    any_fail_rec(r:end) = true;
    break;
end

if early_break_all
    Penalty = cons.pen.bigM * cons.pen.lambda.fail_bigM;
    out = struct();
    out.ratios = struct('spring_tau',0,'spring_slender',0,'stroke',0, ...
                        'force_cap',0,'dp_quant',0,'thermal_dT',0, ...
                        'cav_frac',0,'qsat_margin',0);
    out.any_fail = true;
    out.dpq_records    = dpq_records;
    out.Fmax_records   = Fmax_records;
    out.stroke_records = stroke_records;
    out.dT_records     = dT_records;
    out.cav_records    = cav_records;
    out.Qp95_records   = Qp95_records;
    out.logs           = logs;
    return
end


        % ---- Yön zarfı (X,Y → max)
        Fmax_records(r)   = max(rec_metrics.Fmax,   [], 'omitnan');
        stroke_records(r) = max(rec_metrics.stroke, [], 'omitnan');
        dpq_records(r)    = max(rec_metrics.dpq,    [], 'omitnan');
        dT_records(r)     = max(rec_metrics.dT,     [], 'omitnan');
        cav_records(r)    = max(rec_metrics.cav95,  [], 'omitnan');
        Qp95_records(r)   = max(rec_metrics.Qp95,   [], 'omitnan');
        any_fail_rec(r)   = any(rec_metrics.fail);
    end


    % ---- Kayıt agregasyonları
    % Fmax ve stroke en-kötü kayıt
    Fmax_all   = max(Fmax_records,   [], 'omitnan');
    stroke_all = max(stroke_records, [], 'omitnan');
    dT_all     = max(dT_records,     [], 'omitnan');
    cav_all    = max(cav_records,    [], 'omitnan');
    Qp95_all   = max(Qp95_records,   [], 'omitnan');

    % Δp quantile: kayıtlar arası 'max' veya 'CVaR'
    switch lower(cons.dp.agg)
        case 'cvar'
            dpq_all = cvar_from_samples(dpq_records(:), cons.alpha_CVaR_cons);
        otherwise
            dpq_all = max(dpq_records, [], 'omitnan');
    end

    agg_vals = [Fmax_all stroke_all dpq_all dT_all cav_all Qp95_all];
    if any(~isfinite(agg_vals))
        names = {'Fmax','stroke','dpq','dT','cav','Qp95'};
        msg = ['NaN/Inf aggregate: ' strjoin(names(~isfinite(agg_vals)), ', ')];
        logs{end+1} = msg; %#ok<AGROW>
        log_msg('error', msg);
        Penalty = cons.pen.bigM;
         out = struct('ratios',struct(),'any_fail',true, ...
            'dpq_records',dpq_records,'Fmax_records',Fmax_records, ...
            'stroke_records',stroke_records,'dT_records',dT_records, ...
            'cav_records',cav_records,'Qp95_records',Qp95_records, ...
            'logs',{logs});
        return;
    end

    % ---- Oranlar (≥1 → ihlal)
    ratios = struct();

   % K1: yay kesme gerilmesi (yalnız yay kolu kuvveti)
if cons.on.spring_tau
    k_p_est = sh.G*sh.d_w^4 / (8*sh.n_turn*sh.D_m^3);           % [N/m]
    Fspring_max = k_p_est * stroke_all;                          % [N]
    ratios.spring_tau = spring_tau_ratio(Fspring_max, sh, cons.spring.tau_allow);
else
    ratios.spring_tau = 0;
end

    % K2: burkulma/serbest boy oranı
    if cons.on.spring_slender
        if cons.spring.use_fixed_length && isfinite(cons.spring.L_free_fixed)
            L_free = cons.spring.L_free_fixed;
        else
            L_free = cons.spring.L_free_auto_fac * geom.Lgap; % yaklaşıklama
        end
        lambda = (L_free / max(sh.D_m,eps)) / max(cons.spring.lambda_max,eps);
        ratios.spring_slender = lambda;
    else
        ratios.spring_slender = 0;
    end

    % K3: strok
    if cons.on.stroke
        ratios.stroke = stroke_all / max(cons.stroke.util_factor*geom.Lgap, eps);
    else
        ratios.stroke = 0;
    end

    % K4: cihaz kuvvet limiti
    if cons.on.force_cap
        ratios.force_cap = Fmax_all / max(cons.force.F_cap, eps);
    else
        ratios.force_cap = 0;
    end

    % K5: Δp quantile
    if cons.on.dp_quant
        ratios.dp_quant = dpq_all / max(num.dP_cap, eps);
    else
        ratios.dp_quant = 0;
    end

    % K6: termal ΔT
    if cons.on.thermal_dT
        ratios.thermal_dT = dT_all / max(cons.thermal.cap_C, eps);
    else
        ratios.thermal_dT = 0;
    end

    % K7: kavitasyon
    if cons.on.cav_frac
        ratios.cav_frac = cav_all / max(cons.hyd.cav_frac_cap, eps);
    else
        ratios.cav_frac = 0;
    end

    % K8: Q satürasyon marjı
    if cons.on.qsat_margin
        Qcap = getfield_default(num,'Qcap_big', ...
    0.4 * ( max(max(orf.CdInf,orf.Cd0)*orf.Ao_eff, 1e-9) * sqrt(2*1.0e9 / max(therm.rho_ref,100)) ) );
ratios.qsat_margin = Qp95_all / max(cons.hyd.Q_margin * Qcap, eps);

    else
        ratios.qsat_margin = 0;
    end

    % ---- Ceza hesabı
    hinge = @(r) max(0, r - 1);
    pwr   = cons.pen.power;

    pen = 0;
    if cons.on.spring_tau,     pen = pen + cons.pen.lambda.spring_tau     * hinge(ratios.spring_tau    )^pwr; end
    if cons.on.spring_slender, pen = pen + cons.pen.lambda.spring_slender * hinge(ratios.spring_slender)^pwr; end
    if cons.on.stroke,         pen = pen + cons.pen.lambda.stroke         * hinge(ratios.stroke        )^pwr; end
    if cons.on.force_cap
                               pen = pen + cons.pen.lambda.force_cap      * hinge(ratios.force_cap     )^pwr; end
    if cons.on.dp_quant,       pen = pen + cons.pen.lambda.dp_quant       * hinge(ratios.dp_quant      )^pwr; end
    if cons.on.thermal_dT,     pen = pen + cons.pen.lambda.thermal_dT     * hinge(ratios.thermal_dT    )^pwr; end
    if cons.on.cav_frac,       pen = pen + cons.pen.lambda.cav_frac       * hinge(ratios.cav_frac      )^pwr; end
    if cons.on.qsat_margin,    pen = pen + cons.pen.lambda.qsat_margin    * hinge(ratios.qsat_margin   )^pwr; end

    any_fail = any(any_fail_rec);
    if cons.on.fail_bigM && any_fail
        pen = pen + cons.pen.bigM * cons.pen.lambda.fail_bigM;
    end

    Penalty = pen;

    % Çıkış detayları
    out = struct();
    out.ratios    = ratios;
    out.any_fail  = any_fail;
    out.dpq_records = dpq_records;
    out.Fmax_records = Fmax_records;
    out.stroke_records = stroke_records;
    out.dT_records = dT_records;
    out.cav_records = cav_records;
    out.Qp95_records = Qp95_records;
    out.logs = logs;
end
function ratio = spring_tau_ratio(Fmax, sh, tau_allow)
    C  = max(sh.D_m, eps) / max(sh.d_w, eps);
    Kw = (4*C - 1)/(4*C - 4) + 0.615/C;
    tau_max = (8 * Fmax * max(sh.D_m,eps) / max(pi*sh.d_w^3, eps)) * Kw;
    ratio = tau_max / max(eps, tau_allow);
end

function P = lhs_population(lb,ub,N)
    D = numel(lb); P = zeros(N,D);
    for d=1:D
        edges = linspace(0,1,N+1);
        centers = (edges(1:end-1)+edges(2:end))/2;
        P(:,d) = lb(d) + centers(randperm(N))' .* (ub(d)-lb(d));
    end
end

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
    tern(obj.use_scaled_for_goal,'scaled','raw'), obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num, ...
    cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg), ...
    design_set, x);

% aux içine koy
aux_J1 = J1_split; aux_J2 = J2_split;

% --- Pareto günlüğüne yaz ---
global PARETO;
if isempty(PARETO) || ~isstruct(PARETO)
    PARETO = struct('J1',[],'J2',[],'F',[],'Pen',[],...
                    'set',[],'x',{{}},'feas',[]);
end
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

function cfg2 = cfg_with_ga(design_set,x,geom,sh,orf,hyd,therm,num,cfg)
    % simulate() tasarımı içerde decode ediyor; burada yalnız cfg’yi aynen geçiriyoruz.
    %#ok<INUSD>
    cfg2 = ensure_cfg_defaults(cfg);  % fills PF.auto_t_on, etc.
end
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

% =============================== Alt Yapı ===============================
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


function [t2,a2] = regrid_to_target(t1,a1,prep)
    % Tekil zaman düzelt, hedef dt'ye göre (auto/off/force)
    [t1,iu]=unique(t1,'stable'); a1=a1(iu);
    dt1 = median(diff(t1),'omitnan');
    tgt = prep.target_dt; tol = max(prep.tol_rel*max(tgt,eps), 1e-12);
    switch lower(prep.resample_mode)
        case 'off'
            t2=t1; a2=a1;
        case 'force'
            t2 = (t1(1):tgt:t1(end)).';
            a2 = interp1(t1,a1,t2,prep.regrid_method,'extrap');
        otherwise % 'auto'
            if abs(dt1 - tgt) <= tol
                t2=t1; a2=a1;
            else
                t2 = (t1(1):tgt:t1(end)).';
                a2 = interp1(t1,a1,t2,prep.regrid_method,'extrap');
                warning('Resample: dt=%.6g→%.6g s | N=%d', dt1, tgt, numel(t2));
            end
    end
end

function [t,a,t5,t95] = pick_series(src, rid, dir, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl)

    src = lower(string(src)); dir = upper(string(dir));
    switch src
        case "raw"
            if dir=="X", t=t_rawX{rid}; a=a_rawX{rid}; t5=t5x_raw(rid); t95=t95x_raw(rid);
            else,         t=t_rawY{rid}; a=a_rawY{rid}; t5=t5y_raw(rid); t95=t95y_raw(rid); end
        otherwise % 'scaled'
            if dir=="X", t=t_sclX{rid}; a=a_sclX{rid}; t5=t5x_scl(rid); t95=t95x_scl(rid);
            else,         t=t_sclY{rid}; a=a_sclY{rid}; t5=t5y_scl(rid); t95=t95y_scl(rid); end
    end
    if isempty(a), error('Kayıt #%d için %s yönü mevcut değil.', rid, dir); end
end

function [M,K,C] = make_KCM(n,mv,kv,cv)
    M = diag(mv); K = zeros(n); C = zeros(n);
    for i=1:n
        kL=kv(i); cL=cv(i);
        if i<n, kU=kv(i+1); cU=cv(i+1); else, kU=0; cU=0; end
        K(i,i) = kL + (i<n)*kU;   C(i,i) = cL + (i<n)*cU;
        if i>1, K(i,i-1)=-kL; C(i,i-1)=-cL; end
        if i<n, K(i,i+1)=-kU; C(i,i+1)=-cU; end
    end
end

function [x,a_rel] = lin_MCK_consistent(t, ag, M, C, K)
    n  = size(M,1); r  = ones(n,1);
    dt = median(diff(t));
    agf = griddedInterpolant(t,ag,'linear','nearest');
    odef = @(tt,z) [ z(n+1:end); M \ ( -C*z(n+1:end) - K*z(1:n) - M*r*agf(tt) ) ];
    z0 = zeros(2*n,1);
    opts = odeset('RelTol',2e-3,'AbsTol',1e-6,'MaxStep',max(dt*10,2e-3),'InitialStep',max(dt*0.25,1e-3));
    sol = ode23tb(odef,[t(1) t(end)],z0,opts);
    t_end = sol.x(end); idx = find(t <= t_end + 1e-12);
    if isempty(idx), x=nan(numel(t),n); a_rel=x; warning('lin_MCK_consistent: early stop'); return; end
    t_use = t(idx); Z = deval(sol,t_use).';
    x_use = Z(:,1:n); v_use = Z(:,n+1:end);
    a_use = ( -(M\(C*v_use.' + K*x_use.')).' - ag(1:numel(t_use)).*r.' );
    x=nan(numel(t),n); a_rel=x; x(1:numel(t_use),:)=x_use; a_rel(1:numel(t_use),:)=a_use;
end


  
function [J1, J2] = compute_objectives_split( ...
    src, obj, tail_sec, ...
    t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
    t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
    M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, design_set, x_ga)

% Compute objective components separately:
%   J1 - IDR zarf oranı CVaR
%   J2 - mutlak ivme (RMS+p95) zarf oranı CVaR
    h_story_m = 3.0 * ones(n-1,1); % zaten üstte de var; buraya da koyduk
    [J1, ~] = compute_J1_IDR_over_records( ...
        src, obj, h_story_m, tail_sec, ...
        t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
        t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
        M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, design_set, x_ga);

    [J2, ~] = compute_J2_ACC_over_records( ...
        src, obj, h_story_m, tail_sec, ...
        t_rawX,t_rawY,a_rawX,a_rawY,t_sclX,t_sclY,a_sclX,a_sclY, ...
        t5x_raw,t95x_raw,t5y_raw,t95y_raw,t5x_scl,t95x_scl,t5y_scl,t95y_scl, ...
        M,Cstr,K,n,geom,sh,orf,hyd,therm,num,cfg, design_set, x_ga);
end

%% ---------------------------------------------------------------------
function key = hash_psa_prep(t_cells, a_cells, params)
%HASH_PSA_PREP Create hash for PSA preparation inputs
    buf = [];
    for i = 1:numel(t_cells)
        if ~isempty(t_cells{i}), buf = [buf; t_cells{i}(:)]; end
        if ~isempty(a_cells{i}), buf = [buf; a_cells{i}(:)]; end
    end
    buf = [buf; params(:)];
    key = md5(num2str(buf(:)', '%.16g,'));
end
