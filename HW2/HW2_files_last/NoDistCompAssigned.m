%Example: disturbance compensation
clear;clc;close all

Bd=[1;0];
Cd=eye(2); %for constraints
Cdy=[1 0]; %output matrix for cost function

%Initial condition
xnow=[1;1];
xnext=xnow;
ynext=Cdy*xnow;
ypred=ynext;
%Horizons 
ny=5;
nu=4;

simhor=8;
u=[];
for k=2:simhor+1
    %Update linearization matrices
    x1=xnow(1);x2=xnow(2);
    Ad=[]; %fill in with the Jacobian at x1, x2
       
%Build prediction matrices for output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=Cdy;D=0;
n=size(Ad,1); %state dimension
m=size(Bd,2); %input dimension
p=size(C,1); %output dimension

%Compute Py
P=[];
P=C*Ad; 
for i=1:ny-1
  P=[P; ];
end
Py=P;
%Compute Hy
H=[];
H=zeros(p*ny,m*nu);
for i=1:ny,
    for j=1:i,
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m)=;
    end
    H()=D;
end
%Retain only the first nu blocks
H=H(:,1: );
Hy=H;
    
%Matrices for unit box constraints on states
C=Cd;D=[0;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pc=size(C,1); %output dimension
p=pc;
%Compute Pc
P=[];
P=C*Ad; 
for i=1:ny-1
  P=[P; ];
end
Pc=P;
%Compute Hc
H=[];
H=zeros(p*ny,m*nu);
for i=1:ny,
    for j=1:i,
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m)=;
    end
    H()=D;
end
%Retain only the first nu blocks
H=H(:,1:nu*m);
Hc=H;

%Weight
lambda=2;

%Output reference
r=zeros(ny,1);

%Form the L vectors for output constraints
Ly=repmat(eye(pc),[ny 1]);

Ybar=[1;1];
y_max=Ly*Ybar;
y_min=-y_max;

%Control input limits
UMAX=0.5; 
UB=UMAX*ones(nu,1);
LB=-UB;

S=Hy'*Hy+lambda*eye(m*nu); %Hessian portion

x=xnow;

f=2*Hy'*(Py*x-r);

dyp=y_max-Pc*x;
dym=-y_min+Pc*x;

M=[Hc;-Hc];
d=[dyp;dym];

%Solve quadratic program
ctrlvec=quadprog(S,f,M,d,[],[],LB,UB);
u_apply=ctrlvec(1:m); %extract first term of optimal sequence
ypred_end=Py*xnow+Hy*ctrlvec;
ypred_end=ypred_end(end);
ypred=[ypred ypred_end];
u=[u u_apply];  %store control history
xnext(:,k)=Ad*xnow+Bd*u(:,end); %update plant
ynext=[ynext Cdy*xnext(:,k)]; %store output history
xnow=xnext(:,k); 
end

plot(xnext(1,:),xnext(2,:))
figure(2)
plot(u)
figure(3)
plot(ynext);hold on;plot(ypred,'r')