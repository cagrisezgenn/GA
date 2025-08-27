function P = lhs_population(lb,ub,N)
    D = numel(lb); P = zeros(N,D);
    for d=1:D
        edges = linspace(0,1,N+1);
        centers = (edges(1:end-1)+edges(2:end))/2;
        P(:,d) = lb(d) + centers(randperm(N))' .* (ub(d)-lb(d));
    end
end
