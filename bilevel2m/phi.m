function p=phi(w,n,m)

y=w(n+1:n+m);
z=w(n+m+1:n+2*m);

p=sum(y+z-sqrt(y.^2+z.^2));

end