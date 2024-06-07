function g=gradfScholtes(w)

x=w(1:2);
y=w(3:4);
g(1:8,1)=0;
g(1,1)=-2*x(1);
g(2,1)=-6*x(2);
g(3,1)=-4;
g(4,1)=2*y(2);

end