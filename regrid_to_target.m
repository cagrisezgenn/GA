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
