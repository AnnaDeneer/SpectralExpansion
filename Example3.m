%% Example 3: Glycolytic Oscillator

%% Start-up

clear; clc; close all;
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');
t_PCEstart = tic;

%% Settings for the model

lb     = 0.1;       % Lower bound of the uniformly distributed variable
ub     = 0.5;       % Upper bound
mu     = 0.3;       % Mean of the lognormally distributed variable
sigma  = 0.1;       % Variance
y0     = [1 1];     % Initial conditions
nSteps = 80;        % Number of timepoints
NVar   = 2;         % Number of variables
time   = linspace(0, 100, nSteps); % Linearly spaced timepoints from 0 to 100
odefun = @glyco_osc_ode;           % Function handle to the ODEs

%% Construct matrix for Legendre and Hermite polynomials

N = 10;     % Expansion order

fprintf('==============\n')
fprintf('PCE order: %i \n', N)

% Legendre
Bl = zeros(N,N);
for n = 1:N
    for p = 1:N
        Bl(n,p) = ((n-1)/(sqrt(2*n-1) * sqrt(2*p-1))) * (n==(p+1)) + ...
                 ((p-1)/(sqrt(2*n-1) * sqrt(2*p-1))) * (n==(p-1));
    end
end

% Hermite
Bh = zeros(N,N);
for n = 0:(N-1)
    for p = 0:(N-1)
        Bh(n+1,p+1) = sqrt(n) * (n==(p+1)) + sqrt(p) * (n==(p-1)); 
    end
end

%% Eigenvalues and Eigenvectors for both matrices

% Legendre
[u_nl, eigl] = eig(Bl);
[u_nl,~] = qr(u_nl); % orthonormalize eigenvectors

% Hermite
[u_nh, eigh] = eig(Bh);
[u_nh,~] = qr(u_nh);

%% Legendre Polynomials, univariate
n = 0:N-1;
d = (2.*n+1);
normconst = sqrt(d);

Pl = {};
for i = 1:N
    Pl{i} = flipud(LegendrePoly(i-1)) * normconst(i);
end

%% Hermite Polynomials, univariate
Ph = {};
for i = 1:N
    Ph{i} = flipud(HermitePoly(i-1)') / sqrt(factorial(i-1));
end

%% Calculate l_n functions for Legendre
l_nl=zeros(N,N); % initialize at l_n = 0;
for n = 1:N
    for p = 1:N
        l_nl(n, 1:p) = l_nl(n, 1:p) + (u_nl(p, n) .* Pl{p})';
    end
end

%% Calculate l_n functions for Hermite
l_nh=zeros(N,N); % initialize at l_n = 0;
for n = 1:N
    for p = 1:N
        l_nh(n, 1:p) = l_nh(n, 1:p) + (u_nh(p, n) .* Ph{p})';
    end
end

%% Tensor product of l_n functions
midx{:} = repmat({1:N}, 1, NVar);
midx    = multicomb(midx{:});    % Get the multi-indices
L_n     = cell(1, size(midx,1)); % Initialize cell array for tensor products
for i = 1:size(midx,1)
    comb = midx(i,:);
    L_n{i} = kron(l_nl(comb(1),:), l_nh(comb(2),:));
end

%% Solve ODEs for eigenvalues
lambda_l = diag(eigl);     % Eigenvalues for the uniform variable
lambda_h = diag(eigh);     % Eigenvalues for the lognormal variable
v{1}{1} = lambda_l;
v{1}{2} = lambda_h;
lambda_n = multicomb(v{:});  % All combinations of eigenvalues for both variables

x_n = zeros(numel(n), numel(time));
al = sqrt(log(1+ sigma.^2./mu.^2));
for i = 1:size(lambda_n,1)
    k(1) = (ub+lb)./2 + (ub-lb)./2 .* lambda_n(i,1); % Isoprobabilistic transform uniform
    k(2) = mu.*exp(al.*lambda_n(i,2) - al.^2./2);    % Isoprobabilistic transform lognormal
    [~, y] = ode45(odefun, time, y0, [], k);
    x_n(i,:) = y(:,end)';
end

%% Product of first component of eigenvectors
u_n1l = u_nl(1,:);
u_n1h = u_nh(1,:);
v{1}{1} = u_n1l;
v{1}{2} = u_n1h;
u_n1 = multicomb(v{:});
u_n1prod = prod(u_n1,2);

%% Construct projection
x_tilde = zeros(numel(time), N^NVar);
for i = 1:numel(time)
    for j = 1:numel(u_n1prod)
        x_tilde(i,:) = x_tilde(i,:) + u_n1prod(j).*x_n(j,i).*L_n{j};
    end
end
tPCEproj = toc(t_PCEstart);
fprintf('PCE projection time elapsed: %.4f sec\n', tPCEproj)

%% Compare PCE to numerical solution
% Different parameter settings for (non)-oscillating behaviour:
% oscillation
X(1) = -1;
X(2) = 1.5;
% damped oscillation
% X(1) = -0.5;
% X(2) = 1.5;
% stable fixed point
% X(1) = 1;
% X(2) = 1.5;

% Evaluate polynomial at time points
pce_ex = zeros(1, numel(time));
for j = 1:numel(time)
    pce_ex(j) = multi_polyeval(X, x_tilde(j,:), midx);
end

% Solve ODEs
k(1) = (ub+lb)./2 + (ub-lb)./2 .* X(1); % Isoprobabilistic transform uniform
k(2) = mu.*exp(al.*X(2) - al.^2./2);    % Isoprobabilistic transform lognormal
[t, y] = ode45(odefun, time, y0, [], k);
y = y(:,end);

% Plot comparison
plot(t,y)
hold on; plot(t, pce_ex, '--')
legend('$Y_{true}$', '$Y_{PCE}$')
ylabel('$Y$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')

%% Sensitivity: Sobol Indices (SI)

k1 = find(midx(:,1) > 1 & midx(:,2) == 1);  % Subset of multi-indices that contain first variable
k2 = find(midx(:,2) > 1 & midx(:,1) == 1);  % Subset " second variable
D1 = zeros(1, nSteps);
D2 = zeros(1, nSteps);

for i = 1:nSteps
    D = sum(x_tilde(i,:).^2);
    D1(i) = sum(x_tilde(i,k1).^2)/D;    % First-order SI for first variable
    D2(i) = sum(x_tilde(i,k2).^2)/D;    % First-order SI for second variable
end

% Plot SI over time
figure; bar(time, [D1(:) D2(:)], 'grouped')
ylabel('$\hat{S}_{\theta_i}$', 'interpreter', 'latex');
xlabel('t')
legend('$\theta_1$', '$\theta_2$');
set(gca, 'fontsize', 12)