function g=gradfex922(w)

x=w(1);
y=w(2);
g(1:8,1)=0;
g(1,1)=2*x;
g(2,1)=2*(y-10);
end