function [F, Q, S] = VPML(X, c, alpha, M, N, K, niter, r)
% X is the data matrix
% niter is the number of the maximun iteration
% c is the number of the neighbors when initializing S^v
% M is the number of the views
% N is the number of the data points
% K is the number of the clusters
% r is used to control the condition of the convergence

if nargin < 7
    niter = 30;
end

if nargin < 8
    r = 2;
end

%% initial beta and S^v
beta = cell(M, 1);
S = cell(M, 1);
for p = 1:M
    distX = L2_distance_1(X{p}, X{p});
    [distX1, idx] = sort(distX, 2);
    S{p} = zeros(N);
    for i = 1:N
        di = distX1(i,2:c+2);
        id = idx(i,2:c+2);
        S{p}(i,id) = (di(c+1)-di)/(c*di(c+1) - sum(di(1:c)) + eps);
    end
    beta{p} = alpha;
end




%% variables
bound = 1e-9;
objArr = zeros(niter+1, 1);

%% iteration start
for iter = 1:niter
%     fprintf('Iteration %d...\n', iter);
    
    %% Upate Q
    Q = cell(M, 1);
    for p = 1:M
        W = (S{p}+S{p}')/2;
        Lap{p} = diag(sum(W))-W;
        Q{p} = my_eig(Lap{p}, K, 0);
    end
    
    %% Update F
    Lap = cell(M, 1);
    for p=1:M
        W = (S{p}+S{p}')/2;
        Lap{p} = diag(sum(W))-W;
    end
    Lw = Laplacian_comb(Lap, M, N);
    F = my_eig(Lw, K, 0);
    
    %% Update S
    distF = L2_distance_1(F', F');
    for p = 1:M
        distQ = L2_distance_1(Q{p}', Q{p}');
        S{p} = zeros(N);
        for i = 1:N
            v = (distQ(i,:) + alpha * distF(i,:))/(2*beta{p});
            S{p}(i, :) = SimplexProj(-v);
        end
    end   
    objArr(iter) = objfun(F, Q, S, M, alpha, beta);


    
    %% Update beta, alpha
    alpha = 0;
    E = cell(M,1);
    for p = 1:M
        E{p} = (S{p}+S{p}')/2;
        D = diag(sum(E{p}));
        L = D-E{p};
        [~, ~, ev]=my_eig(L, c, 0);
        fn1 = sum(ev(1:c));
        fn2 = sum(ev(1:c+1));
        if fn1 > 0.00000000001
            beta{p} = beta{p}*0.5;
        elseif fn2 < 0.00000000001
            beta{p} = beta{p}/0.5;
        else
            break;
        end
        alpha = alpha + beta{p};
    end
    alpha = alpha/M;
    
    %% stopping criterion
    if iter > 1
        if abs((objArr(iter)-objArr(iter-1))/(N^r)) < bound
            break;
        end
    end
end

end

