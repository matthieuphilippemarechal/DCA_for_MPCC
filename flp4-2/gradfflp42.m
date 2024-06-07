function g=gradfflp42(w,n,m)

x=w(1:n,1);
g(n+2*m,1)=0;
g(1:n,1)=x;
g(n+1:n+m,1)=ones(1,m);

end