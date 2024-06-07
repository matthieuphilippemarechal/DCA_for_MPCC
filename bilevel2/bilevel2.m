%%%%%% Solves min f(w) with w=(x,y,z,l) Aeq*w=beq, A*w<=b, l>=0, z>=0, z.l=0
clear all
p=8;
m=12;
M=30;

A(1,1:32)=0;
A(1,1:4)=[1,1,1,1];
b=40;

Aeq(1:16,1:32)=0;
Aeq(1,5)=1;
Aeq(1,9:12)=[0.4,0.6,-1,1];
Aeq(2,6)=1;
Aeq(2,9:14)=[0.7,0.3,0,0,-1,1];
Aeq(3,7)=1;
Aeq(3,15:18)=[0.4,0.6,-1,1];
Aeq(4,8)=1;
Aeq(4,15:20)=[0.7,0.3,0,0,-1,1];


Aeq(5,1)=1;
Aeq(5,5:6)=[-0.4,-0.7];
Aeq(5,21)=-1;
Aeq(6,2)=1;
Aeq(6,5:6)=[-0.6,-0.3];
Aeq(6,22)=-1;
Aeq(7,5)=1;
Aeq(7,23)=-1;
Aeq(8,5)=-1;
Aeq(8,24)=-1;
Aeq(9,6)=1;
Aeq(9,25)=-1;
Aeq(10,6)=-1;
Aeq(10,26)=-1;
Aeq(11,3)=1;
Aeq(11,7)=-0.4;
Aeq(11,8)=-0.7;
Aeq(11,27)=-1;
Aeq(12,4)=1;
Aeq(12,7)=-0.6;
Aeq(12,8)=-0.3;
Aeq(12,28)=-1;
Aeq(13,7)=1;
Aeq(13,29)=-1;
Aeq(14,7)=-1;
Aeq(14,30)=-1;
Aeq(15,8)=1;
Aeq(15,31)=-1;
Aeq(16,8)=-1;
Aeq(16,32)=-1;

beq(1:16,1)=0;
beq(1:4,1)=[4;13;35;2];
beq(8,1)=-20;
beq(10,1)=-20;
beq(14,1)=-40;
beq(16,1)=-40;

borneinf(1:32,1)=0;
bornesup(1:32,1)=100;
bornesup(1:4,1)=[10;5;15;20];
bornesup(4:7,1)=[20;20;40;40];
bornesup(23:26,1)=20;
bornesup(27:28,1)=[15;20];
bornesup(29:32,1)=40;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=0.1;
tolc=0.001;
itermax=1000;
errc=2*tolc;
iterc=0;
itermaxc=1000;
rho=10;

wt(1:32,1)=0;
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
    g2=-gradfbilevel2(w(:,iter))+lambda*w(:,iter);
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
   while j<=q
       wstar=gradPhitilde(w(:,iter+1),p,m,j,I,plop);
%        if w(2)^2+w(3)^2<0.001
%            wstar=[0;0;2];
%        end
       wt=quadprog(eye(p+2*m),gradfbilevel2(w(:,iter+1))-w(:,iter+1)+mu*wstar,A,b,Aeq,beq,borneinf,bornesup);
       errc=norm(wt-w(:,iter+1),Inf);       
       if errc>tolc
           j=q+1;
       else
           j=j+1;
       end
       rho=fbilevel2(wt)-fbilevel2(w(:,iter+1))+mu*(phi(wt,p,m)-phi(w(:,iter+1),p,m));
   end
       ww(:,iterc)=w(:,iter+1);
       iiter(iterc)=iter+1;       
end

% ww(:,iterc+1)=wt;
 fval=fbilevel2(wt);
 
 x=wt(1:4);
 y=wt(5:8);
 l=wt(9:20);
 z=wt(21:32);
