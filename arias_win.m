function [t5,t95] = arias_win(t,ag,p1,p2)
    if nargin<3 || isempty(p1), p1=0.05; end
    if nargin<4 || isempty(p2), p2=0.95; end
    t=t(:); ag=ag(:);
    mask=isfinite(t) & isfinite(ag);
    t=t(mask); ag=ag(mask);
    if numel(t)<2, t5=t(1); t95=t(end); return; end
    [t,ia]=unique(t,'stable'); ag=ag(ia);
    dt=diff(t); a2=ag.^2;
    Eincr=0.5*(a2(1:end-1)+a2(2:end)).*dt;
    E=[0; cumsum(Eincr)]; Eend=E(end);
    if Eend<=eps, t5=t(1); t95=t(end); return; end
    t5=taf(p1); t95=taf(p2);
    t5=max(min(t5,t(end)),t(1));
    t95=max(min(t95,t(end)),t(1));
    if t95<=t5, t5=t(1); t95=t(end); end
    function tv=taf(frac)
        target=frac*Eend;
        k=find(E>=target,1,'first');
        if isempty(k)
            tv=t(end);
        elseif k==1
            tv=t(1);
        else
            E0=E(k-1); E1=E(k);
            r=(target-E0)/max(E1-E0,eps);
            r=min(max(r,0),1);
            tv=t(k-1)+r*(t(k)-t(k-1));
        end
    end
end
