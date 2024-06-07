function xj=basedeux(j,plop)
xj=0*ones(1,plop);
while j~=0
   r=mod(j,2);
   xj(plop)=r;
   plop=plop-1;
   j=(j-r)/2;
end

end