function dydt = glyco_osc_ode(t, y, k)
dydt=zeros(size(y));
x = y(1);
Y = y(2);
a = k(1);
b = k(2);
dydt(1) = -x + a*Y + x^2*Y;
dydt(2) = b - a*Y - x^2*Y;

end