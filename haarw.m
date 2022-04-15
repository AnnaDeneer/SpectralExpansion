function psi = haarw(x, j, k)
psi = zeros(1, numel(x));

for i = 1:numel(x)
    if j < 0
        psi(i) = 1;
    elseif (2^j * x(i) - k) >= 0 && (2^j * x(i) - k) < 1/2
        psi(i) = 2^(j/2)*1;
    elseif (2^j * x(i) - k) > 1/2 && (2^j * x(i) - k) <= 1
        psi(i) = 2^(j/2)*-1;
    else 
        psi(i) = 0;
    end
end

end