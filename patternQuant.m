function TD = patternQuant(k)

TS = lsa(k);
if TS
    % initialization
    NVar = 3;
    ctr = IJKth(1,1:20,1:20,20,NVar);
    [~, y] = runModel(k);
    ssAC = y(end, ctr+2);
    T = find(ssAC >= 0.5*max(ssAC));
    nT = numel(T);
    TD = nT/400;  % trichome density
else
    TD = 0;
end

end

