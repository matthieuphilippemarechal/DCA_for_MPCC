function g=fbilevel2m(w)


y=w(5:8);
g=- (200 - y(1) - y(3))*(y(1) + y(3)) - (160 - y(2) - y(4))*(y(2) + y(4));
end