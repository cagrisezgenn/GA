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
