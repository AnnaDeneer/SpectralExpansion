function dxdt = odefunc(t,x,k)
    dxdt = zeros(size(x));  
    
    dxdt(1) = k(1) - k(2)*x(1)*x(2) + k(3)*x(3) - k(5)*x(1);
    dxdt(2) = k(4) - k(2)*x(1)*x(2) + k(3)*x(3) - k(5)*x(2);
    dxdt(3) = k(2)*x(1)*x(2) - k(3)*x(3) - k(5)*x(3);
end