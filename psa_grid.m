function [tPSA, agPSA] = psa_grid(t, ag, down_dt)
%PSA_GRID Optional resampling for PSA (downsample only; no upsample).
%   [tPSA, agPSA] = PSA_GRID(t, ag, down_dt) resamples the acceleration
%   record AG defined at times T onto a coarser grid with spacing DOWN_DT.
%   If DOWN_DT is empty or not greater than the native time step, the input
%   vectors are returned unchanged. Resampling uses PCHIP interpolation and
%   extrapolates outside the original domain.

    t = t(:); ag = ag(:);
    if nargin < 3 || isempty(down_dt) || down_dt <= 0
        tPSA = t; agPSA = ag; return;
    end
    dt0 = median(diff(t), 'omitnan');
    if down_dt <= dt0 * (1 + 1e-12)
        tPSA = t; agPSA = ag; return;
    end
    tPSA  = (t(1):down_dt:t(end)).';
    agPSA = interp1(t, ag, tPSA, 'pchip', 'extrap');
end
