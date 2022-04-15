function [t, y] = runModel(k, tspan, varargin)
% Script to run the trichome patterning model

if nargin < 2
    tspan = logspace(0, 3, 100);
end

xmax = 20; % number of cells in x-direction
ymax = 20; % " in y-direction

% Depletion model
model = @SD;
NVar = 3;
% label = {'TTG1', 'GL3', 'AC'};
%% Initial conditions
% Initialize cell grid indices
ctr = IJKth(1,1:ymax,1:xmax,ymax,NVar);
% Construct spatial coupling matrix for hexagonal cells and periodic
% boundary conditions, last arg: 1 = zero flux, 0 = periodic
D = couplingMatrix(ymax,xmax,[-1 1 0 0 1 -1],[0 0 -1 1 -1 1],1);

% Set integration options including structure of the jacobian
options = odeset('Vectorized','on','JPattern',jpat(NVar,D));

% Compute homogeneous (=single cell) steady state: start from zero
% protein levels and leave sufficient time to reach steady state
% Use 1 % random perturbations of the homogeneous steady state as
% initial conditions for the tissue model
% [~,y] = ode15s(model, tspan, ones(NVar,1), [], 1, 0, k);
% ss =y(end, :);
% y0 = repmat(ss(:),ymax*xmax,1) .* (1 + 0.01.*rand(NVar*ymax*xmax,1));
y0 = load('y0trichome.mat'); y0 = y0.y0;

%% Integrate the model
[t,y] = ode15s(model, tspan, y0(:), options, ctr, D, k);



end