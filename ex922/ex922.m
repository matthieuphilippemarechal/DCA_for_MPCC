%%%%%% Solves min f(x,y,z) with A[x;y]=z, B[x;y]<=c, y>=0, z>=0, y.z=0
clear all
n=2;
m=3;
M=10;

A(1,1:8)=0;

A(1,1:2)=[-1,1];

b=0;

Aeq(1:4,1:8)=0;
Aeq(1,1:3)=[1,1,1];
Aeq(2,1:4)=[0,-1,0,1];
Aeq(3,1:5)=[0,1,0,0,1];
Aeq(4,1:2)=[2,4];
Aeq(4,6:8)=[1,-1,1];


beq=[20;0;20;60];

borneinf(1:8,1)=0;
bornesup(1:8,1)=100;
bornesup(1:2,1)=15;
bornesup(3,1)=20;
bornesup(4,1)=15;
bornesup(5,1)=20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=0.1;
tolc=0.0001;
itermax=1000;
err=2*tol;
err2=2*tolc;
iter=0;
errc=2*tolc;
iterc=0;
itermaxc=100;
rho=10;

wt(1:8,1)=3;
aalpha=10;
mu=170;
lambda=20;
while (errc>tolc)&&(iterc<itermaxc)
    clear w
    iterc=iterc+1;
    w(:,1)=wt;
    iter=0;
while (err>tol)&&(iter<itermax)
    iter=iter+1;    
    g2=-gradfex922(w(:,iter))+lambda*w(:,iter);
g2=g2-mu*gradPhi(w(:,iter),n,m,rho/(M*(2^iter)));
[w(:,iter+1),f,exitflag,output,mult]=quadprog(lambda*eye(n+2*m),-g2,A,b,Aeq,beq,borneinf,bornesup);
lambdamult=mult.eqlin;
mumult=mult.ineqlin;
betainf=mult.lower;
betasup=mult.upper;
err=norm(gradfex922(w(:,iter+1))+mu*gradPhi(w(:,iter),n,m,rho/(M*(2^iter)))+A'*mumult+Aeq'*lambdamult-betainf+betasup);
%err=norm(w(:,iter+1)-w(:,iter));
end
y=w(n+1:n+m,iter+1);
z=w(n+m+1:n+2*m,iter+1);
plop=0;
I=[];
for j=1:m
    if y(j)^2+z(j)^2<=0.1  
     plop=plop+1;
     I=[I,j];
    end
end
q=2^(plop);
   j=1;
   while j<=q
       wstar=gradPhitilde(w(:,iter+1),n,m,j,I,plop);
%        if w(2)^2+w(3)^2<0.001
%            wstar=[0;0;2];
%        end
       wt=quadprog(eye(n+2*m),gradfex922(w(:,iter+1))-w(:,iter+1)+mu*wstar,A,b,Aeq,beq,borneinf,bornesup);
       errc=norm(wt-w(:,iter+1));       
       if errc>tolc
           j=q+1;
       else
           j=j+1;
       end
       rho=fex922(wt)-fex922(w(:,iter+1))+mu*(phi(wt,n,m)-phi(w(:,iter+1),n,m));
   end
       ww(:,iterc)=w(:,iter+1);
       iiter(iterc)=iter+1;
end

 ww(:,iterc+1)=wt;
 fval=fex922(wt);
 
 x=wt(1);
 y=wt(2);
 s=wt(3:5);
 l=wt(6:8);
