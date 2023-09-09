% code for 2013 arXiv paper
% Projection onto the probability simplex: An efficient algorithm with a simple proof, and an application
function X = SimplexProj(Y)
	[N,D] = size(Y);
	X = sort(Y,2,'descend');
	Xtmp = (cumsum(X,2)-1)*diag(sparse(1./(1:D)));
	X = max(bsxfun(@minus,Y,Xtmp(sub2ind([N,D],(1:N)', sum(X>Xtmp,2)))),0);
end