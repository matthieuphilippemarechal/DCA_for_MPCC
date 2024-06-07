%%%%%% Solves min f(x,y,z) with A[x;y]=z, B[x;y]<=c, y>=0, z>=0, y.z=0
clear all
n=50;
m=30;
Ms=300;

[AA,b]=matrixA(1);
M=matrixM(1);
N=matrixN(1);
beq=-matrixq(1);

A=[AA,0*ones(30,60)];

Aeq=[N,M,-eye(30)];

borneinf(1:n,1)=-Inf;
borneinf(n+1:n+2*m,1)=0;
bornesup(1:n+2*m,1)=1000;

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
itermaxc=1000;
rho=10;

wt(1:n+2*m,1)=0;
aalpha=10;
mu=1000;
lambda=300;
while (errc>tolc)&&(iterc<itermaxc)
    clear w
    iterc=iterc+1;
    w(:,1)=wt;
    iter=0;
while (err>tol)&&(iter<itermax)
    iter=iter+1;    
    g2=-gradfflp42(w(:,iter),n,m)+lambda*w(:,iter);
g2=g2-mu*gradPhi(w(:,iter),n,m,rho/(Ms*(2^iter)));
w(:,iter+1)=quadprog(lambda*eye(n+2*m),-g2,A,b,Aeq,beq,borneinf,bornesup);
err=norm(w(:,iter+1)-w(:,iter));
end
y=w(n+1:n+m,iter+1);
z=w(n+m+1:n+2*m,iter+1);
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
       wstar=gradPhitilde(w(:,iter+1),n,m,j,I,plop);
        wt=quadprog(eye(n+2*m),-w(:,iter+1)+gradfflp42(w(:,iter+1),n,m)+mu*wstar,A,b,Aeq,beq,borneinf,bornesup);
%       wt=linprog(gradfBard1(w(:,iter+1))+mu*wstar,[],[],Aeq,beq,borneinf,bornesup);
%       errc=(gradfBard1(w(:,iter+1))+mu*wstar)'*(w(:,iter+1)-wt);
       errc=norm(wt-w(:,iter+1));       
       if errc>tolc
           j=q+1;
       else
           j=j+1;
       end
       rho=fflp42(wt,n,m)-fflp42(w(:,iter+1),n,m)+mu*(phi(wt,n,m)-phi(w(:,iter+1),n,m));
   end
       ww(:,iterc)=w(:,iter+1);
end

 %ww(:,iterc+1)=wt;
 fval=fflp42(wt,n,m);
 
 x=wt(1:n);
 y=wt(n+1:n+m);
 z=wt(n+m+1:n+2*m);
