%% Example 5: Trichome model, direct expansion

%% Start-up
clear; clc; close all;
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');


%% Model settings
% Parameter values and standard deviation taken from Bouyer et al. 2008
% PLoS Biology

k       = [0.6662 0.1767 3.1804 5.3583]; % Default parameter values
pidx    = 1;                             % Index of stochastic variable
sigma   = [0.2690 0.0677 1.4742 2.4240]; % Standard deviation of parameters
lb      = k(pidx) - sigma(pidx);         % Lower bound of uniform distribution
ub      = k(pidx) + sigma(pidx);         % Upper bound
NVar    = numel(pidx);                   % Number of stochastic variables


%% Construct matrix
N = 8; % Expansion order
fprintf('==============\n')
fprintf('PCE order: %i \n', N)
B = zeros(N,N);

% Legendre
for n = 1:N
    for p = 1:N
        B(n,p) = ((n-1)/(sqrt(2*n-1) * sqrt(2*p-1))) * (n==(p+1)) + ...
                 ((p-1)/(sqrt(2*n-1) * sqrt(2*p-1))) * (n==(p-1));
    end
end
%% eigenvalues and eigenvectors
[u_n, eigvals] = eig(B);
[u_n,r] = qr(u_n); % orthonormalize eigenvectors


%% Legendre Polynomials, univariate
n = 0:N-1;
d = (2.*n+1);
normconst = sqrt(d);

P = {};
for i = 1:N
    P{i} = flipud(LegendrePoly(i-1)) * normconst(i);
end


%% Calculate l_n functions
t_PCEstart = tic;
l_n=zeros(N,N); % initialize at l_n = 0;
for n = 1:N
    for p = 1:N
        l_n(n, 1:p) = l_n(n, 1:p) + (u_n(p, n) .* P{p})';
    end
end

%% Evaluate model at eigenvalues
lambda_n = diag(eigvals);
x_n = load('x_nDirect.mat'); x_n = x_n.x_n;

% Uncomment to calculate what is loaded from .mat file:
% x_n = zeros(N, 1);
% 
% for i = 1:N
%     k(pidx) = (ub+lb)./2 + (ub-lb)./2 .* lambda_n(i); % Isoprobabilistic transform (uni)
%     x_n(i) = patternQuant(k);
% end

%% Construct projection
x_tilde = zeros(1, N^NVar);

for j = 1:N
    x_tilde = x_tilde + u_n(1,j).*x_n(j).*l_n(j,:);
    
end

tPCEtot = toc(t_PCEstart);
fprintf('PCE construction time elapsed: %.4f sec\n', tPCEtot)

%% Compare PCE to numerical solutions

X = linspace(-1, 1, 10^2)'; % Range of samples
% Ytrue = zeros(1,size(X,1));
for j = 1:size(X,1)
    k = [0.6662 0.1767 3.1804 5.3583];
    x = X(j);
    p =  (ub+lb)./2 + (ub-lb)./2 .* x; % Isoprobablistic transform
    k(pidx) = p;
    if lsa(k)
        Ypce(j)  = polyval(fliplr(x_tilde), x);
%         Ytrue(j) =  patternQuant(k);
    else
        Ypce(j) = 0;
%         Ytrue(j) = 0;
    end
end
% Ytrue above is commented-out and loaded from .mat file for speed
Ytrue = load('Ytrue_trichome.mat'); Ytrue = Ytrue.Ytrue;

% Kernel density estimation
[f_mc, xi_mc]   = ksdensity(Ytrue);
[f_pce, xi_pce] = ksdensity(Ypce);

figure;
plot(xi_mc, f_mc, 'k'); hold on;
plot(xi_pce, f_pce, '--');
legend('$Y_{true}$', '$Y_{PCE}$')
ylabel('Probability density', 'interpreter', 'latex')
xlabel('Y', 'interpreter', 'latex')


