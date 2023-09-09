function [Y, J] = AWSC(Q, S, truthY, niter, VN)
% VN-th view is selected to initalize Y. Generally, the view which has the
% largest dimension is suggested to be selected.
M = length(S);
N = size(S{1}, 1);
K = length(unique(truthY));

%% Initalize Y
CKSym = (S{VN}+S{VN}')/2;
predY = SpectralClustering(CKSym, K);
Y = zeros(K, N);
for i = 1:N
    Y(predY(i),i) = 1;
end

%% Initalize p_v
for v = 1:M
    p(v) = 1/M;
    Q{v} = Q{v}';
end

G = cell(M, 1);
J = zeros(niter+1, 1);
for iter = 1:niter
%     fprintf('Iteration %d...\n', iter);
    
    %% Update G^v
    for v = 1:M
        G{v} = Q{v}*transpose(Y)/(Y*transpose(Y));
    end
    
    %% Update Y
    Y = zeros(K, N);
    E = eye(K);
    for i = 1:N
        Sum = zeros(K, 1);
        for j = 1:K
            for v = 1:M
                Sum(j) = Sum(j) + norm(Q{v}(:,i)-G{v}*E(:,j), 2)^2/p(v);
            end
        end
        k = find(Sum == min(Sum));
        Y(k,i) = 1;
    end
    
    %% Update p_v
    for v = 1:M
        p(v) = norm(Q{v}-G{v}*Y, 'fro');
    end
    p = p/sum(p);
    %% stopping 
    J(iter) = Jfun(Q, Y, G, p);
%     fprintf('objective = %d\n', J(iter));
    if iter > 1
        if abs(J(iter)-J(iter-1)) < 1e-12
            break;
        end
    end
end

end

