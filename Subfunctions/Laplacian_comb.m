function Lw = Laplacian_comb(L, M, N)
w = zeros(M,1);
Lw = zeros(N, N);
maxiter = 20;
eps = 1e-12;
%initalize weight
for p = 1:M
    Lw = Lw + L{p}/M;
end
for iter = 1:maxiter
    for p = 1:M
        w(p) = 1/(2*norm(Lw-L{p},'fro'));
    end
    w = w./sum(w);
    for p = 1:M
        Lw = Lw + w(p)*L{p};
    end
    obj(iter) = 0;
    for p = 1:M
        obj(iter) = obj(iter)+norm(Lw-L{p},'fro');
    end
    if iter > 1
        if abs(obj(iter)-obj(iter-1)) < eps
            break;
        end
    end
end
end

