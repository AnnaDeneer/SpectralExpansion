function dydt = schnak(t, y, ctr, D, k)

u = y(ctr,:);
v = y(ctr+1,:);
a = k(1);
b = k(2);   
gamma = k(3);
d = k(4);

dydt(ctr,:) = gamma*(a - u + u.^2.*v) + (D*u);

dydt(ctr+1,:) = gamma*(b - u.^2.*v) + d.*(D*v); 

end