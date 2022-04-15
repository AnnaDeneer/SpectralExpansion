%% Example 5: Trichome model, indirect expansion

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
Ncells  = 400;                           % Number of cells on the grid

%% construct matrix
N = 8;
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

%% Evaluate model at lambda_n
lambda_n = diag(eigvals);
x_n = load('x_nIndirect.mat'); x_n = x_n.x_n;

% Uncomment to calculate what is loaded from .mat file:
% x_n = zeros(N, 400);
for i = 1:N
    k(pidx) = (ub+lb)./2 + (ub-lb)./2 .* lambda_n(i); % Isoprobabilistic transform (uni)
    x_n(i,:) = patternQuant_ind(k);
end

%% Construct projection
x_tilde = zeros(Ncells, N^NVar);
for i = 1:Ncells
    for j = 1:N
        x_tilde(i,:) = x_tilde(i,:) + u_n(1,j).*x_n(j,i).*l_n(j,:);
    end
end
tPCEtot = toc(t_PCEstart);
fprintf('PCE construction time elapsed: %.4f sec\n', tPCEtot)

%% Comparison of PCE to numerical solution
X = linspace(-1, 1, 10^3)';
Ypce   = zeros(size(X,1), Ncells);
Ytrue  = zeros(size(X,1), Ncells);
TDpce  = zeros(size(X,1), 1);
TDtrue = zeros(size(X,1) , 1);

% Evaluate PCE
pces = tic;
for j = 1:size(X,1)
    k = [0.6662 0.1767 3.1804 5.3583];
    x = X(j);
    p =  (ub+lb)./2 + (ub-lb)./2 .* x; 
    k(pidx) = p;
    for c = 1:Ncells
        Ypce(j,c) = polyval(fliplr(x_tilde(c,:)), x);
    end
    % Post-processing of PCE evaluation:
    if lsa(k) % Check for Turing Space
        % The trichome density is determined by finding those points
        % where the concentration is higher than half-maximum
        TDpce(j) = numel(find(Ypce(j,:) >= 0.5*max(Ypce(j,:))))/400;
    else
        TDpce(j) = 0;
    end
end
tpce = toc(pces);
fprintf('Time elapsed %i samples PCE: %.2f sec\n', size(X,1), tpce)

% Evaluate true model
% truet = tic;
% for j = 1:size(X,1)
%     k = [0.6662 0.1767 3.1804 5.3583];
%     x = X(j);
%     p =  (ub+lb)./2 + (ub-lb)./2 .* x; 
%     k(pidx) = p;
%     Ytrue(j,:) = patternQuant_ind(k);
%     % Post-processing of numerical solution:
%     if lsa(k)
%         TDtrue(j) = numel(find(Ytrue(j,:) >= 0.5*max(Ytrue(j,:))))/400;
%     else
%         TDtrue(j) = 0;
%     end
% end
% ttrue = toc(truet);
% fprintf('Time elapsed %i samples true model: %.2f sec\n', size(X,1), ttrue)

% Uncomment above to calculate true solutions given in the .mat here:
TDtrue = load('TDtrue_indirect.mat'); TDtrue = TDtrue.TDtrue;

% Kernel density estimation
[f_mc, xi_mc] = ksdensity(TDtrue);
[f_pce, xi_pce] = ksdensity(TDpce);

figure;
plot(xi_mc, f_mc, 'k'); hold on;
plot(xi_pce, f_pce, '--');
legend('$Y_{true}$', '$Y_{PCE}$')
ylabel('Probability density', 'interpreter', 'latex')
xlabel('Y', 'interpreter', 'latex')


