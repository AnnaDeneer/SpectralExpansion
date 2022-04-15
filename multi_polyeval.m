function y = multi_polyeval(x, coefs, terms)
pow = terms-1;
NVar = numel(x);
y = 0;
for i = 1:size(pow,1)
    yc = coefs(i); % coefficient
    for j = 1:NVar
        yc = yc * x(j)^pow(i,j); % x_j^p
    end
    y = y + yc;
end




end