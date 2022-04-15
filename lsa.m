function TS = lsa(k)
%% Determine single cell stability

TS = false;
NVar = 3;
xmax = 20;
ymax = 20;
dpq = FourierModes(xmax,ymax, 'hex');
% ctr = IJKth(1,1:xmax,1:ymax,ymax,NVar);
model = @SD;
[~,y] = ode15s(model, [0 10000], ones(NVar,1), [], 1, 0, k);
ss = y(end,:);
TTG1 = ss(1); GL3 = ss(2); AC = ss(3);
J = [ - GL3 - k(2),      -TTG1,  0;...
       -GL3, - TTG1 - 1, 2*AC*k(4);...
        GL3,       TTG1,      -1];
D = diag([k(3), 0, 0]);
%Find positive eigenvalues for J:
pos_eig = find(eig(J)>0);
if (numel(pos_eig) >= 1)
    return
end

%% Check stability for multi-cell Apq matrix

for idx = length(dpq):-1:1
    Apq = J-D*dpq(idx);
    % Find positive eigenvalues:
    evals = eig(Apq);
    pos_eig = find(evals>0);
    
    % If for any given p,q the matrix Apq has an eigenvalue
    % with a positive real part the ss is unstable:
    if (numel(pos_eig) >= 1)
        TS = true;
        break;
    end
end

% if TS
%     [~, y] = runModel(k, logspace(0, 2, 100));
%     ssAC = y(end, ctr+2);
%     if var(ssAC) < (10^-6)
%         TS = false;
%     end
% end
