function Jvalue = Jfun(Q, Y, G, p)
% The function of AWSC
M = length(p);
Jvalue = 0;
for v = 1:M
    Jvalue = Jvalue + norm(Q{v}-G{v}*Y, 'fro')^2/p(v);
end

end

