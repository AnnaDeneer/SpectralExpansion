    function dydt = SD(t, y, ctr, D, k)%, z)
    % Trichome patterning model extension of activator inhibitor & depletion
    % including dimers (Bouyer et al. 2008 PLoS Biology)
    dydt = zeros(size(y));

    TTG1 = y(ctr,:);
    GL3 = y(ctr+1,:);
    AC = y(ctr+2,:); % TTG1 - GL3

    % TTG1
    dydt(ctr,:)   = k(1) - TTG1.*k(2) - TTG1.*GL3...
                    + k(3).*(D*TTG1);

    % GL3
    dydt(ctr+1,:) = k(4).*AC.^2 - GL3 - TTG1.*GL3;

    % AC: TTG1 - GL3
    dydt(ctr+2,:) = TTG1.*GL3 - AC ;
    end