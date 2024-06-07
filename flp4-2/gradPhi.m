function Y=gradPhi(w,n,m,rho)

Y(1:n+2*m)=0;

y=w(n+1:n+m);
z=w(n+m+1:n+2*m);

for k=1:m
   if y(k)+z(k)+rho>0.001
       Y(n+k)=1-y(k)/(sqrt(y(k)^2+z(k)^2+rho^2));
       Y(n+m+k)=1-z(k)/(sqrt(y(k)^2+z(k)^2+rho^2));
   else
       r=rand(1);
%       r=1; 
       theta=(pi/2)*rand(1);
       Y(n+k)=1-r*cos(theta);
       Y(n+m+k)=1-r*sin(theta);
   end

end
Y=Y';
end