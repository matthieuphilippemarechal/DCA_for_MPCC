%%%%%% Solves min f(w) with w=(x,y,z,l) Aeq*w=beq, A*w<=b, l>=0, z>=0, z.l=0
clear all
p=4;
m=2;
M=30;

Aeq(1:4,1:8)=0;
Aeq(1,3:6)=[2,0,2,-3];
Aeq(2,5:6)=[-1,4];
Aeq(3,1:4)=[-2,0,-2,1];
Aeq(3,7)=-1;
Aeq(4,2:4)=[1,3,-4];
Aeq(4,8)=-1;

beq(1:4,1)=[0;5;-3;4];

A(1,1:8)=0;
A(1,2)=2;

b=4;

borneinf(1:8,1)=0;
bornesup(1:8,1)=100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=0.1;
tolc=0.001;
itermax=1000;
errc=2*tolc;
iterc=0;
itermaxc=100;
rho=10;
aalpha=1000;
mu=100;
lambda=300;

wt(1:8,1)=0;

while (errc>tolc)&&(iterc<itermaxc)
    clear w
    iterc=iterc+1;
    w(:,1)=wt;
    iter=0;
    err=2*tol;
while (err>tol)&&(iter<itermax)
    iter=iter+1;    
    A(1,1)=2*w(1,iter);
    Aeq(3,1)=2*w(1,iter)-2;
    Aeq(3,2)=2*w(2,iter);
    b=4+w(1,iter)^2;
    beq(3)=-3+w(1,iter)^2-w(2,iter)^2;
    g2=-gradfBard3(w(:,iter))+lambda*w(:,iter);
g2=g2-mu*gradPhi(w(:,iter),p,m,rho/(M*(2^iter)));
w(:,iter+1)=quadprog(lambda*eye(p+2*m),-g2,A,b,Aeq,beq,borneinf,bornesup);
err=norm(w(:,iter+1)-w(:,iter),Inf);
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
       A(1,1)=2*w(1,iter+1);
    Aeq(3,1)=2*w(1,iter+1)'-2;
    Aeq(4,1)=2*w(2,iter+1);
    b=4+w(1,iter+1)^2;
    beq(3)=-3+w(1,iter+1)^2+w(2,iter+1)^2;
   while j<=q
       wstar=gradPhitilde(w(:,iter+1),p,m,j,I,plop);
        wt=quadprog(eye(p+2*m),-w(:,iter+1)+gradfBard3(w(:,iter+1))+mu*wstar,A,b,Aeq,beq,borneinf,bornesup);
%       wt=linprog(gradfBard3(w(:,iter+1))+mu*wstar,[],[],Aeq,beq,borneinf,bornesup);
%       errc=(gradfBard3(w(:,iter+1))+mu*wstar)'*(w(:,iter+1)-wt);
       errc=norm(wt-w(:,iter+1),Inf);       
       if errc>tolc
           j=q+1;
       else
           j=j+1;
       end
       rho=fBard3(wt)-fBard3(w(:,iter+1))+mu*(phi(wt,p,m)-phi(w(:,iter+1),p,m));
   end
       ww(:,iterc)=w(:,iter+1);
end

 %ww(:,iterc+1)=wt;
 fval=fBard3(wt);
 
 x=wt(1:2);
 y=wt(3:4);
 l=wt(5:6);
 z=wt(7:8);
