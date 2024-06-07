%%%%%% Solves min f(w) with w=(x,y,z,l) Aeq*w=beq, A*w<=b, l>=0, z>=0, z.l=0
clear all
p=2;    
m=3;
M=30;

A(1,1:8)=0;

A(1,1:5)=-[-1.5,2,1,-0.5,1];
b=-2;

Aeq(1:3,1:8)=0;
Aeq(1,1:2)=[3,-1];
Aeq(1,6)=-1;
Aeq(2,1:2)=[-1,0.5];
Aeq(2,7)=-1;
Aeq(3,1:2)=[-1,-1];
Aeq(3,8)=-1;

beq(1:3,1)=0;
beq(1:3,1)=[3;-4;-7];

borneinf(1:8,1)=0;
bornesup(1:8,1)=100;
bornesup(1:2,1)=7;
bornesup(6:8,1)=[18;7.5;7];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=0.001;
tolc=0.001;
itermax=1000;
err=2*tol;
err2=2*tolc;
iter=0;
errc=2*tolc;
iterc=0;
itermaxc=1000;
rho=10;

wt(1:8,1)=0;
aalpha=10;
mu=1000;
lambda=30;
while (errc>tolc)&&(iterc<itermaxc)
    clear w
    iterc=iterc+1;
    w(:,1)=wt;
    iter=0;
    err=2*tol;
while (err>tol)&&(iter<itermax)
    iter=iter+1;    
    g2=-gradfBard1m(w(:,iter))+lambda*w(:,iter);
g2=g2-mu*gradPhi(w(:,iter),p,m,rho/(M*(2^iter)));
w(:,iter+1)=quadprog(lambda*eye(p+2*m),-g2,A,b,Aeq,beq,borneinf,bornesup);
err=norm(w(:,iter+1)-w(:,iter));
end
y=w(p+1:p+m,iter+1);
z=w(p+m+1:p+2*m,iter+1);
plop=0;
I=[];
for j=1:m
    if y(j)^2+z(j)^2<=0.01  
     plop=plop+1;
     I=[I,j];
    end
end
q=2^(plop);
   j=1;
   while j<=q
       wstar=gradPhitilde(w(:,iter+1),p,m,j,I,plop);
%        if w(2)^2+w(3)^2<0.001
%            wstar=[0;0;2];
%        end
       wt=linprog(gradfBard1m(w(:,iter+1))+mu*wstar,A,b,Aeq,beq,borneinf,bornesup);
       errc=norm(wt-w(:,iter+1));       
       if errc>tolc
           j=q+1;
       else
           j=j+1;
       end
       rho=fBard1m(wt)-fBard1m(w(:,iter+1))+mu*(phi(wt,p,m)-phi(w(:,iter+1),p,m));
   end
       ww(:,iterc)=w(:,iter+1);
end

 ww(:,iterc+1)=wt;
 fval=fBard1m(wt);
 
 x=wt(1);
 y=wt(2);
 l=wt(3:5);
 z=wt(6:8);
