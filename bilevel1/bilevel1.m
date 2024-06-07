%%%%%% Solves min f(w) with w=(x,y,z,l) Aeq*w=beq, A*w<=b, l>=0, z>=0, z.l=0
clear all
p=4;
m=6;
M=10;
A(1,1:16)=0;
A(1,1:4)=[1,1,1,-2];
Aeq(1:8,1:16)=0;
Aeq(1:2,1:10)=[-2,0,2,0,-1,1,0,0,2,0;0,-2,0,2,0,0,-1,1,0,2];
Aeq(3,3)=1;
Aeq(3,11)=-1;
Aeq(4,3)=-1;
Aeq(4,12)=-1;
Aeq(5,4)=1;
Aeq(5,13)=-1;
Aeq(6,4)=-1;
Aeq(6,14)=-1;
Aeq(7,1)=1;
Aeq(7,3)=-2;
Aeq(7,15)=-1;
Aeq(8,2)=1;
Aeq(8,4)=-2;
Aeq(8,16)=-1;

b=0;
beq=[-40;-40;-10;-20;-10;-20;10;10];
d(1:16,1)=0;
d(1:4,1)=[2;2;-3;-3];

borneinf(1:16,1)=0;
borneinf(3:4,1)=-10;
bornesup(1:2,1)=50;
bornesup(3:4,1)=20;
bornesup(5:16,1)=inf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=0.1;
tolc=0.001;
itermax=1000;
wt(1:16,1)=0;
wt(1:4,1)=[0;0;0;0];
err=2*tol;
err2=2*tolc;
iter=0;
aalpha=10;
mu=50;
lambda=10;
errc=2/tolc;
iterc=0;
itermaxc=100;
rho=10;
while (errc>tolc)&&(iterc<itermaxc)
    clear w
    iterc=iterc+1;
    w(:,1)=wt;
    iter=0;
    err=2*tol;
while (err>tol)&&(iter<itermax)
    iter=iter+1;    
    g2=-d+lambda*w(:,iter);
g2=g2-mu*gradPhi(w(:,iter),p,m,rho/(M*(2^iter)));
w(:,iter+1)=quadprog(lambda*eye(p+2*m),-g2,A,b,Aeq,beq,borneinf,bornesup);

err=norm(w(:,iter+1)-w(:,iter),Inf);
end
y=w(p+1:p+m,iter+1);
z=w(p+m+1:p+2*m,iter+1);
plop=0;
I=[];
for j=1:m
    if y(j)^2+z(j)^2<=0.001  
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
       wt=quadprog(eye(p+2*m),d-w(:,iter+1)+mu*wstar,A,b,Aeq,beq,borneinf,bornesup);
       errc=norm(wt-w(:,iter+1),Inf);       
       if errc>tolc
           j=q+1;
       else
           j=j+1;
       end
       rho=d'*wt-d'*w(:,iter+1)+mu*(phi(wt,p,m)-phi(w(:,iter+1),p,m));
   end
       ww(:,iterc)=w(:,iter+1);
       iiter(iterc)=iter+1;
end

 ww(:,iterc+1)=wt;
 fval=d'*wt-60;
 
 x=wt(1:2);
 y=wt(3:4);
 l=wt(5:10);
 z=wt(11:16);
