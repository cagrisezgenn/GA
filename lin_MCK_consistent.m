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
