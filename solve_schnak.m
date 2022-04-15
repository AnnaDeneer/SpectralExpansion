function [t,y] = solve_schnak(k)
model = @schnak;
NVar = 2;
N = 20; % number of cells
D = couplingMatrix(1,N,[-1 1 0 0],[0 0 -1 1],1); %square grid, zero flux
ctr = IJKth(1,1,1:N,1,NVar);
tspan = linspace(0, 10^6, 1000);

a = k(1); b = k(2);
u0 = a + b;
v0 = b/(a+b)^2;
ss = [u0 v0];

% Initial conditions
ry0 = load('y0schnak.mat'); ry0 = ry0.ry0;
y0 = repmat(ss(:),N,1) .* (1 + 0.01.*ry0);

options = odeset('Vectorized','on','JPattern',jpat(NVar,D));
[t,y] = ode15s(model, tspan, y0(:), options, ctr, D, k);

end