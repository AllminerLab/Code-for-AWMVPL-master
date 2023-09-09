clear
close all


addpath('Subfunctions')
%% load data
%% Mfeat
load('Mfeat.mat');
dataname = 'Mfeat';
alpha = 0.5;
c = 25;
VN = 5;
type = 1;
r = 3;

fprintf('Datasets: %s, c = %d, alpha = %f ...\n', dataname, c, alpha);

%% parameters
K = length(unique(y));
M = length(X);
N = length(y);
niter = 30;

%% Normalization
X = Normalization(X,type);

%% VPML
[F, Q, S] = VPML(X, c, alpha, M, N, K, niter, r);

%% AWSC 
[Y, J] = AWSC(Q, S, y, niter, VN);


%% Show the clustering performance
preY = zeros(N, 1);
for i = 1:N
    preY(i) = find(Y(:, i) == 1);
end
[acc,nmi,pur] = ClusteringMeasure(y, preY);
fprintf('Clusteing performance:\n');
fprintf('ACC = %.3f, NMI = %.3f, PUR = %.3f\n', acc, nmi, pur);
fprintf('------------------------------------------\n');
fprintf('\n');

