%% Example 4: Schnakenberg model

%% Start-up
clear; clc; close all;
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

%% Model settings
lb = 0.001;         % Lower bound of uniform distribution
ub = 0.45;          % Upper bound
NVar = numel(lb);   % Number of stochastic variables
NSpecies = 2;       % Number of species
Ncells   = 20;      % Number of cells on the grid
ctr      = IJKth(1,1,1:Ncells,1,NSpecies); % Discritization indices

%% Concatenate wavelet scaling and sliding factors
N = 6;                     % Wavelet resolution level
x = linspace(-1, 1, 1000); % Haar wavelets are defined over [-1, 1]
F = @(x) 1/2*(1+x);        % Transform
Nw = (2^N-1)*2 + 1;        % Total nr of wavelets for N res. lvls
Nwa = 128;                 % Sub-selection
psi = zeros(Nw, numel(x));
idx = zeros(Nw, 2);

for j = 0:N % Scaling factor
    % Uncomment for wavelet plots:
%     subplot(1,N+1,j+1);
    for k = 0:2^j-1 % Sliding factor
        idx(2^j+k,:) = [j k];
        psi(2^j+k,:) = haarw(F(x),j,k);
        % Uncomment for wavelet plots:
%         plot(x,psi(2^j+k,:), 'linewidth', 2); hold on;
%         title(['j = ' num2str(j)])
    end
end
idx = [-1 0; idx];

%% Calculate matrix for Haar wavelets
fprintf('==============\n')
fprintf('Wavelet resolution: %i \n', N)
P = 1/2; % Weight function of uniform distribution
A = zeros(Nwa,Nwa);
for n = 1:Nwa
    for m = 1:Nwa
        f = @(x) haarw(F(x),idx(n,1),idx(n,2)) .* haarw(F(x),idx(m,1),idx(m,2)).*P .* x;
        A(n,m) = integral(f, -1, 1);
    end
end

%% Eigenvalues and eigenvectors
[u_n, eigvals] = eig(A);
[u_n,r] = qr(u_n); % orthonormalize eigenvectors

%% Evaluate model at eigenvalues
tic_lambda = tic;
lambda_n = diag(eigvals);
x_n = zeros(size(lambda_n,1), Ncells);

fprintf('Solving model for %i eigenvalues...\t', numel(lambda_n))
for i = 1:Nwa
    k = (ub+lb)./2 + (ub-lb)./2 .* lambda_n(i); % Isoprobabilistic transform
    pars = [k 1 5 20];     % Fix other parameters
    [~,y] = solve_schnak(pars);
    x_n(i,:) = y(end,ctr); % Get steady state solution for all cells
end

tlambda = toc(tic_lambda);
fprintf('Time elapsed: %.4f sec\n', tlambda)

%% Comparison of PCE v. real model solution
x = 0.006;                            % Chosen parameter
k = (ub+lb)./2 + (((ub-lb)./2) .* x); % Transform
pars = [k 1 5 20];                    % Fix other parameters
Ypce = zeros(1, Ncells);

tPCEsum = tic;
% Summation
for j = 1:Ncells
    for l = 1:Nwa
        for n = 1:Nwa
            Ypce(j) = Ypce(j) + x_n(l,j) * u_n(1,l) * u_n(n,l) * (haarw(F(x), idx(n,1), idx(n,2)));
        end
    end
end
tPCE = toc(tPCEsum);
fprintf('PCE summation time elapsed: %.4f sec\n', tPCE)

% Solve discretized ODEs
tODEstart = tic;
[~,y] = solve_schnak(pars);
ytrue = y(end,ctr);
tODE = toc(tODEstart);
fprintf('ODE solving time elapsed: %.4f sec\n', tODE)

% Plot comparison
figure;
plot(1:Ncells,ytrue, 'k', 'linewidth', 1)
hold on; plot(1:Ncells, Ypce, '--.', 'Color', [0.64 0.38 0.38],...
              'markersize', 10, 'linewidth', 1)
legend(['$Y_{true}, \alpha =' num2str(k, '%.2f') '$'], '$Y_{Haar}$');
ylabel('$Y$', 'interpreter', 'latex')
xlabel('$cell$', 'interpreter', 'latex')


