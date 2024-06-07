function g=fflp42(w,n,m)

x=w(1:n,1);
y=w(n+1:n+m,1);
g=0.5*x'*x+ones(1,m)*y;
end