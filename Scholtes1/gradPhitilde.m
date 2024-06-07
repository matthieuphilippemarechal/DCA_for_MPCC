function Y=gradPhitilde(w,n,m,j,I,plop)

xj=basedeux(j,plop+1);

Y(1:n+2*m)=0;

y=w(n+1:n+m);
z=w(n+m+1:n+2*m);

for k=1:m
   if y(k)+z(k)>0.1
       Y(n+k)=1-y(k)/(sqrt(y(k)^2+z(k)^2));
       Y(n+m+k)=1-z(k)/(sqrt(y(k)^2+z(k)^2));
   end
end

for k=1:plop
    if xj(k)==0
        Y(n+I(k))=1;
    else
        Y(n+m+I(k))=1;
    end

end
Y=Y';
end