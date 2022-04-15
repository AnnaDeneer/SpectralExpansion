function dpq = FourierModes(xmax, ymax, gridtype)
dpq = zeros(1,xmax*ymax);
idx = 0;
for p = 1:xmax
    for q = 1:ymax
        idx = idx +1;
        if strcmp(gridtype, 'square')
             dpq(idx) = sin(pi*p/xmax)^2 + sin(pi*q/ymax)^2;
        else %hexagrid
        dpq(idx) = 4*(sin(pi*p/xmax)^2 + sin(pi*q/ymax)^2 + ...
            sin((pi*p/xmax)-(pi*q/ymax))^2);
        end
    end
end
dpq = unique(dpq);
end