function obj = objfun(F, Q, S, M, alpha, beta)
% The function of VPML
obj = 0;
    for p = 1:M
        W = (S{p}+S{p}')/2;
        L = diag(sum(W)) - W;
        obj = obj + (2*trace(Q{p}'*L*Q{p}) + 2*alpha*trace(F'*L*F) + beta{p} * norm(S{p}, 'fro').^2);
    end
end

