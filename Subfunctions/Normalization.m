function X = Normalization(X,type)
% The input X \in R^{n*d^{v}}, the output X \in R^{d^{v}*n}
% If X is a non-sparse data matrix, type is suggected to be 1
% If X is a sparse data matrix, type is suggected to be 2
M = length(X);
N = size(X{1},1);
if type == 1
    for p = 1:M
        X{p} = X{p}'./ norm(X{p}, 'fro');
    end
elseif type == 2
        for v = 1:M
            for i = 1:N
                normItem = std(X{v}(i,:));
                if (0 == normItem)
                    normItem = eps;
                end
                X{v}(i,:) = (X{v}(i,:) - mean(X{v}(i,:)))/normItem;
            end
            X{v} = X{v}';
        end
end
end

