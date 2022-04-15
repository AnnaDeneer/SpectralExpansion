%% Example 2: Biochemical Reaction Network

%% Start-up
clear; clc; close all;
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

%% Model settings
nSteps = 20;                        % Number of time points
lb = [0.1 0.1 0.1 0.1 0.1];         % Lower bound of uniform distributions
ub = [1 1 1 1 1];                   % Upper bound of uniform distributions
time = linspace(0, 200, nSteps);    % Linearly spaced time points from 0 to 200
y0 = [0 0 0];                       % Initial conditions
odefun = @binding;                  % Function handle to ODEs
NVar = numel(lb);                   % Number of variables

%% Construct matrix
N = 3; % Expansion order
fprintf('==============\n')
fprintf('PCE order: %i \n', N)
B = zeros(N,N);
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
l_n=zeros(N,N); % initialize at l_n = 0;
for n = 1:N
    for p = 1:N
        l_n(n, 1:p) = l_n(n, 1:p) + (u_n(p, n) .* P{p})';
    end
end

%% Evaluate model at eigenvalues
% Note: can take considerable computation time

tlambda_start = tic;

% Get all combinations of eigenvalues:
v{:} = repmat({diag(eigvals)},1,NVar);
lambda_n = multicomb(v{:});
x_n1 = zeros(numel(n), numel(time));
x_n2 = zeros(numel(n), numel(time));
x_n3 = zeros(numel(n), numel(time));

% Segmentation settings:
M = 1; % Segmentation granularity
L = 1; % Interval length
G = @(k,m) (2*m*L)/(2*M+1) + k/(2*M+1); % Scaling function
m = -M:1:M; % Segment indices
mm{:} = repmat({m},1,NVar);
mm_n = multicomb(mm{:}); % Segment multi-indices

% Get numerical solutions for all segments:
for ii = 1:size(mm_n, 1)
    g = G(lambda_n, mm_n(ii,:));
    for i = 1:size(g,1)
        k = (ub+lb)./2 + (ub-lb)./2 .* g(i,:);
        [~, y] = ode45(odefun, time, y0, [], k);
        x_n1(i,:) = y(:,1)';   % First variable
        x_n2(i,:) = y(:,2)';   % Second variable
        x_n3(i,:) = y(:,end)'; % Last variable
    end
    x_m1{ii} = x_n1;
    x_m2{ii} = x_n2;
    x_m3{ii} = x_n3;
end

tlambda = toc(tlambda_start);
fprintf('Model evaluation time: %.4f sec\n', tlambda)

%% Construct projection
tsumstart = tic;
v{:} = repmat({u_n(1,:)}, 1, NVar);
u_n1 = multicomb(v{:});  % All combinations of first component of eigenvectors
u_n1prod = prod(u_n1,2); % Product of first component

% Tensor product of l_n functions
midx{:} = repmat({1:N}, 1, NVar);
midx = multicomb(midx{:});  % All multi-indices
L_n = cell(1, size(midx,1));
for i = 1:size(midx,1)
    comb = midx(i,:);
    L_n{i} = kron(l_n(comb(1),:), l_n(comb(2),:));
    for j = 3:numel(comb)
        L_n{i} = kron(L_n{i}, l_n(comb(j),:));  
    end
end

% Projection for each variable and each sub-segment:
for ii = 1:size(mm_n, 1)
    x_n1 = x_m1{ii};
    x_n2 = x_m2{ii};
    x_n3 = x_m3{ii};
    x_tilde1 = zeros(numel(time), N^NVar);
    x_tilde2 = zeros(numel(time), N^NVar);
    x_tilde3 = zeros(numel(time), N^NVar);
    for i = 1:numel(time)
        for j = 1:numel(u_n1prod)
            x_tilde1(i,:) = x_tilde1(i,:) + u_n1prod(j).*x_n1(j,i).*L_n{j};
            x_tilde2(i,:) = x_tilde2(i,:) + u_n1prod(j).*x_n2(j,i).*L_n{j};
            x_tilde3(i,:) = x_tilde3(i,:) + u_n1prod(j).*x_n3(j,i).*L_n{j};
        end
    end
    x_tilde_m1{ii} = x_tilde1;
    x_tilde_m2{ii} = x_tilde2;
    x_tilde_m3{ii} = x_tilde3;
end
tsum = toc(tsumstart);
fprintf('Summation time elapsed: %.4f sec\n', tsum)

%% Compare PCE to numerical solution

k = [0.556 0.111 0.9 0.111 0.111]; % Parameter values

% Select corresponding segment:
Z    = @(t) floor(((2*M+1)* t)/(2*L) + 0.5); % Index function
Ginv = @(k,m) (2*M+1) * k - 2*L*m;           % Scaling function
zi   = Z(k);                                 % Indices for chosen parameter values
z    = 1 + (M+zi(1))*(2*M+1)^4 + (M+zi(2))*(2*M+1)^3 + (M+zi(3))*(2*M+1)^2 + ...
       (M + zi(4))*(2*M+1) + (M+zi(5));      % The segment multi-index

% Evaluate polynomials
x_tilde1 = x_tilde_m1{z};
x_tilde2 = x_tilde_m2{z};
x_tilde3 = x_tilde_m3{z};
pce_ex1 = zeros(1, nSteps);
pce_ex2 = zeros(size(pce_ex1));
pce_ex3 = zeros(size(pce_ex1));
for j = 1:numel(time)
        pce_ex1(j) = multi_polyeval(Ginv(k,Z(k)), x_tilde1(j,:), midx);
        pce_ex2(j) = multi_polyeval(Ginv(k,Z(k)), x_tilde2(j,:), midx);
        pce_ex3(j) = multi_polyeval(Ginv(k,Z(k)), x_tilde3(j,:), midx);
end

% Numerical solution using ode45
k = (ub+lb)./2 + (((ub-lb)./2) .* k);
[t, y] = ode45(odefun, time, y0, [], k);

% Plot comparison
figure; plot(t,y(:,1), 'k')
hold on; plot(t, pce_ex1, 'k.')
plot(t,y(:,2), 'r')
plot(t, pce_ex2, 'r.')
plot(t,y(:,3), 'b')
plot(t, pce_ex3, 'b.')
legend('$Y_{true}$', '$Y_{PCE}$')
ylabel('$Y$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')

