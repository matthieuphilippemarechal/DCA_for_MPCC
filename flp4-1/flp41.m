%%%%%% Solves min f(w) with w=(x,y,z) Aeq*w=beq, A*w<=b, y>=0, z>=0, z.y=0
clear all
p=50;
m=30;
Ms=6*m;
[AA,b]=matrixA(1);
M=matrixM(1);
N=matrixN(1);
beq=-matrixq(1);

A=[AA,0*ones(30,60)];

Aeq=[N,M,-eye(30)];

borneinf(1:p,1)=-Inf;
borneinf(p+1:p+2*m,1)=0;
bornesup(1:p+2*m,1)=Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=0.1;
tolc=0.001;
itermax=1000;
err=2*tol;
err2=2*tolc;
iter=0;
errc=2*tolc;
iterc=0;
itermaxc=10000;
rho=10;

wt(1:p+2*m,1)=0;
aalpha=10;
mu=1000;
lambda=1;
while (errc>tolc)&&(iterc<itermaxc)
    clear w
    err=2*tol;
    iterc=iterc+1;
    w(:,1)=wt;
    iter=0;
while (err>tol)&&(iter<itermax)
    iter=iter+1;    
    g2=-gradfflp41(w(:,iter),p,m)+lambda*w(:,iter);
g2=g2-mu*gradPhi(w(:,iter),p,m,rho/(Ms*(2^iter)));
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
   while j<=q
       wstar=gradPhitilde(w(:,iter+1),p,m,j,I,plop);
        wt=quadprog(eye(p+2*m),-w(:,iter+1)+gradfflp41(w(:,iter+1),p,m)+mu*wstar,A,b,Aeq,beq,borneinf,bornesup);
%       wt=linprog(gradfBard1(w(:,iter+1))+mu*wstar,[],[],Aeq,beq,borneinf,bornesup);
%       errc=(gradfBard1(w(:,iter+1))+mu*wstar)'*(w(:,iter+1)-wt);
       errc=norm(wt-w(:,iter+1),Inf);       
       if errc>tolc
           j=q+1;
       else
           j=j+1;
       end
       rho=fflp41(wt,p,m)-fflp41(w(:,iter+1),p,m)+mu*(phi(wt,p,m)-phi(w(:,iter+1),p,m));
   end
       ww(:,iterc)=w(:,iter+1);
end

 %ww(:,iterc+1)=wt;
 fval=fflp41(wt,p,m);
 
 x=wt(1:p);
 y=wt(p+1:p+m);
 z=wt(p+m+1:p+2*m);