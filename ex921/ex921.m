%%%%%% Solves min f(w) with w=(x,y,s,l) Aeq*w=beq, A*w<=b, l>=0, s>=0, s.l=0
clear all
p=2;
m=4;
M=10;

Aeq(1:5,1:10)=0;
Aeq(1,1:3)=[-3,1,1];
Aeq(2,1:4)=[1,-0.5,0,1];
Aeq(3,1:5)=[1,1,0,0,1];
Aeq(4,1:6)=[0,-1,0,0,0,1];
Aeq(5,1:2)=[-1.5,2];
Aeq(5,7:10)=[1,-0.5,1,-1];


beq=[-3;4;7;0;2];

borneinf(1:10,1)=0;
bornesup(1:10,1)=100;
bornesup(1:2,1)=7;
bornesup(3,1)=18;
bornesup(4,1)=8;
bornesup(5,1)=7;
bornesup(6,1)=7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=0.1;
tolc=0.001;
itermax=1000;
err2=2*tolc;
errc=2*tolc;
iterc=0;
itermaxc=100;
rho=10;

wt(1:10,1)=3;
aalpha=10;
mu=170;
lambda=20;
while (errc>tolc)&&(iterc<itermaxc)
    clear w
    iterc=iterc+1;
    w(:,1)=wt;
    iter=0;
    err=2*tol;
while (err>tol)&&(iter<itermax)
    iter=iter+1;    
    g2=-gradfex921(w(:,iter))+lambda*w(:,iter);
g2=g2-mu*gradPhi(w(:,iter),p,m,rho/(M*(2^iter)));
w(:,iter+1)=quadprog(lambda*eye(p+2*m),-g2,[],[],Aeq,beq,borneinf,bornesup);
err=norm(w(:,iter+1)-w(:,iter),Inf);
end
y=w(p+1:p+m,iter+1);
z=w(p+m+1:p+2*m,iter+1);
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
       wstar=gradPhitilde(w(:,iter+1),p,m,j,I,plop);
%        if w(2)^2+w(3)^2<0.001
%            wstar=[0;0;2];
%        end
       wt=quadprog(eye(p+2*m),gradfex921(w(:,iter+1))-w(:,iter+1)+mu*wstar,[],[],Aeq,beq,borneinf,bornesup);
       errc=norm(wt-w(:,iter+1),Inf);       
       if errc>tolc
           j=q+1;
       else
           j=j+1;
       end
       rho=fex921(wt)-fex921(w(:,iter+1))+mu*(phi(wt,p,m)-phi(w(:,iter+1),p,m));
   end
       ww(:,iterc)=w(:,iter+1);
end

 ww(:,iterc+1)=wt;
 fval=fex921(wt);
 
 x=wt(1);
 y=wt(2);
 s=wt(3:6);
 l=wt(7:10);
