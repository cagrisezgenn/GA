function ratio = spring_tau_ratio(Fmax, sh, tau_allow)
    C  = max(sh.D_m, eps) / max(sh.d_w, eps);
    Kw = (4*C - 1)/(4*C - 4) + 0.615/C;
    tau_max = (8 * Fmax * max(sh.D_m,eps) / max(pi*sh.d_w^3, eps)) * Kw;
    ratio = tau_max / max(eps, tau_allow);
end
