function Sab = sdof_PSA_band_avg_aug(t, ag, T1, zeta, band_fac, Np)
% Band ortalaması alınmış spektral ivme
    Tvec = linspace(band_fac(1)*T1, band_fac(2)*T1, Np);
    Sa   = sdof_PSA_vec_aug_ode(t, ag, Tvec, zeta);
    Sab  = mean(Sa);
end
