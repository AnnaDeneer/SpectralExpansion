%% Example 1: Exponential decay

%% Initialization
clear; clc; close all;
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');


%% Settings for the model
nSteps = 13;                     % Number of time points
time   = linspace(0, 6, nSteps); % Linearly spaced time points
y0     = 10;                     % Initial condition
mu     = 0.5;                    % Mean of the lognormal distribution
sigma  = 0.2;                    % Variance
pd = makedist('Lognormal','mu',log(mu),'sigma',sigma);
x = (0:0.01:1.2);
y = pdf(pd, x);

%% Plot the PDF
figure;
z = zeros(size(x));
col = x;  % This is the color, vary with x in this case.
surface([x; x],[y;y],[z;z],[col;col],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
xlabel('k(\theta)')
ylabel('Probability density')
set(gca, 'fontsize', 12)

%% Construct B matrix
N = 5; % Expansion order
fprintf('==============\n')
fprintf('PCE order: %i \n', N)
B = zeros(N,N);
for n = 0:(N-1)
    for p = 0:(N-1)
        B(n+1,p+1) = sqrt(n) * (n==(p+1)) + sqrt(p) * (n==(p-1)); 
    end
end

%% Get eigenvalues and eigenvectors from B matrix
[u_n, eigvals] = eig(B);
[u_n,r] = qr(u_n); % orthonormalize eigenvectors

%% Hermite Polynomials, univariate
P = {};
for i = 1:N
    P{i} = flipud(HermitePoly(i-1)') / sqrt(factorial(i-1));
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
x_n = zeros(N, numel(time));

for i = 1:N
    k = mu.*exp(sqrt(log(1+(sigma^2/mu^2))).*lambda_n(i,:) - ...
        0.5*log(1 + sigma^2/mu^2) ); % Isoprobabilistic transform
    y = y0.*exp(-k.*time);
    x_n(i,:) = y;
end

%% Construct projection
x_tilde = zeros(numel(time), N);
for i = 1:numel(time)
    for j = 1:N
        x_tilde(i,:) = x_tilde(i,:) + u_n(1,j).*x_n(j,i).*l_n(j,:);
    end
end
tPCEtot = toc(t_PCEstart);
fprintf('PCE construction time elapsed: %.4f sec\n', tPCEtot)
%% MC vs PCE
Xval = randn([10^3 1]);  % Generate 1000 MC samples

tidx = 11;               % Choose a single timepoint
t    = time(tidx);         

% Sample PCE for PDF samples
YPCE = zeros(size(Xval, 1), 1);
for i = 1:size(Xval,1)
    k = Xval(i);
    YPCE(i) = polyval(fliplr(x_tilde(tidx,:)), k);
end

% MC for PDF
YM = zeros(size(YPCE));
K = mu.*exp(sqrt(log(1+(sigma^2/mu^2))).*Xval - ...
        0.5*log(1 + sigma^2/mu^2) ); % Isoprobablistic transform
for i = 1:size(Xval,1)
    k = K(i);
    y = y0.*exp(-k.*t);
    YM(i) = y;
end

% Sum of squared error between MC and PCE
err = (YPCE - YM).^2;
serr = sum(err,2);
fprintf('Total squared error: %.4f \n', sum(serr))

% Plot the resulting PDFs
[f_mc, xi_mc,umc]    = ksdensity(YM);   % PDF by MC
[f_pce, xi_pce,upce] = ksdensity(YPCE); % PDF by PCE

figure;
plot(xi_mc, f_mc); hold on;
plot(xi_pce, f_pce, '.'); 

legend('$Y_{true}$', '$Y_{PCE}$')
ylabel('Probability density', 'interpreter', 'latex')
xlabel('Y', 'interpreter', 'latex')
title(['PDF at t = ', num2str(t)])
set(gca, 'fontsize', 12)


