function ssAC = patternQuant_ind(k)

    % initialization
    NVar = 3;
    ctr = IJKth(1,1:20,1:20,20,NVar);
    [~, y] = runModel(k);
    ssAC = y(end, ctr+2);

end

