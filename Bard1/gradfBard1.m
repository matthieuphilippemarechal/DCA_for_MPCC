function g=gradfBard1(w)

x=w(1);
y=w(2);
g(1:8,1)=0;
g(1,1)=2*(x - 5);
g(2,1)=4*(2*y + 1);

end