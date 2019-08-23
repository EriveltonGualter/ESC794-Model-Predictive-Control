%Example: Dual Decomposition Distributed MPC -
%With input constraints and terminal equilibrium constraints

%Programmed by H Richter for MPC course: Selected Topics on Engineering
%Science

%Cleveland State University, Mechanical Engineering Department
%Fall 2018

function ddMPCexampleUconst

clear;clc;close all

% Decoupled plant
Ad = [1 2 0;0 -1 0;0 0 -2];
Bd = [1 0;-1 0;0 1];

%Objective: regulation to an equilibrium state
%Find some eq state for a specific equilibrium control
ueq = [0;0];
xeq = [0; 0; 0];

n = size(Ad,1); %state dimension
p = n;
m = size(Bd,2); %input dimension

% Horizon
ny = 5;
nu = ny;

%Solve and simulate centralized problem
Q = eye(3);
R = eye(2);

[Px, Hx] = buildPH(Ad,Bd,nu,ny,n);

%Set up matrices for constraints on control

Cconstr=eye(3); %dummy
pc=size(Cconstr,1); %dummy
[Pc,Hc,Cc]=buildPHC(Cconstr,Ad,Bd,ny,nu,pc);

xeq_hat=xeq;
ueq_hat=ueq;
Qbar=Q;
Rbar=R;

u_max=[15;15];
u_min=[-15;-15];


ymax=[100;100;100]; %dummy
ymin=[-100;-100;-100]; %dummy
[M,d,Qbar,Rbar,xeq_hat,ueq_hat]=buildMdEq(Cc,xeq,ueq,Q,R,ymax,ymin,u_max,u_min,ny,nu,pc);

%Form cost function - constant portion

S=2*(Hx'*Qbar*Hx+Rbar);
simhor=20; 

%Simulation (unit step reference)
%Initial conditions
xa=[1;1;1];

%initialize variable part of cost
f=2*((xa'*Px'-xeq_hat')*Qbar*Hx-ueq_hat'*Rbar);

uprev=zeros(m,1);
xnow=zeros(3,1);
u=uprev;

%Set up equality constr matrices
Qc=[zeros(n,n*(ny-1)) eye(n)]*Px;
Aeq=[zeros(n,n*(ny-1)) eye(n)]*Hx;
%beq defined inside mpc loop

T=[];
for k=1:simhor;
   beq=xeq-Qc*xa;
   t=tic;
   ctrlvec=quadprog(S,f,M,d,Aeq,beq);
   T(k)=toc(t);
   u_apply=ctrlvec(1:m); %extract first term of optimal sequence
   u=[u u_apply];  %store control history
   xnext(:,k)=Ad*xnow+Bd*u(:,end); %update plant
   xnow=xnext(:,k); 
   uprev=u(:,end);
   xa=xnow;
   f=2*((xa'*Px'-xeq_hat')*Qbar*Hx-ueq_hat'*Rbar);
end

figure(1);
subplot(3,1,1)
stairs(xnext(1,:));hold on;
subplot(3,1,2)
stairs(xnext(2,:));hold on;
subplot(3,1,3)
stairs(xnext(3,:));hold on;


figure(2);
subplot(2,1,1);
stairs(u(1,:));hold on;
subplot(2,1,2)
stairs(u(2,:));hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now solve a dual decomposition MPC with the partition {x1 x2}, {x3}


%Reset plant state
xa=[1;1;1];

%Prepare sub-problems

A11=Ad(1:2,1:2); A22=Ad(3,3); %decoupled subsystem matrices
A12=Ad(1:2,3);A21=Ad(3,1:2); %coupling matrices 

B1=Bd(1:2,1); B2=Bd(3,2); 
Q1c=Qc(:,1:2);Q2c=Qc(:,3);

R1c=Aeq(:,[1 3 5 7 9]); %odd cols of Aeq
R2c=Aeq(:,[2 4 6 8 10]); %even cols of Aeq

%Build prediction matrices for subsystem 1: (A11, B1)

p1=2;
[Px1,Hx1]=buildPH(A11,B1,nu,ny,p1);

%Set up matrices for constraints on control

%Compute Pc, Hc and Cc for constraints

Cd_constr1=eye(2); %dummy
pc1=size(Cd_constr1,1); %dummy

[Pc1,Hc1,Cc1]=buildPHC(Cd_constr1,A11,B1,ny,nu,pc1);
%Pc, Hc unused here

%Compute M and d
xeq1=xeq(1:2);
ueq1=ueq(1);
umax1=u_max(1);
umin1=u_min(1);

R1=1;
Q1=eye(2); %for local optimization


ymax1=[100;100]; %dummy
ymin1=[-100;-100]; %dummy
[M1,d1,Qbar1,Rbar1,xeq1_hat,ueq1_hat]=buildMdEq(Cc1,xeq1,ueq1,Q1,R1,ymax1,ymin1,umax1,umin1,ny,nu,pc1);

%Form cost function - quadratic term 
S1=2*Hx1'*Qbar1*Hx1+Rbar1;




xa1=xa(1:2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Build matrices for subsystem 2: (A12, B2)

p2=1;
[Px2,Hx2]=buildPH(A22,B2,nu,ny,p2);

%Set up matrices for constraints on control

%Compute Pc, Hc and Cc for constraints

Cd_constr2=eye(1); %dummy
pc2=size(Cd_constr2,1); %dummy


[Pc2,Hc2,Cc2]=buildPHC(Cd_constr2,A22,B2,ny,nu,pc2);
%Pc and Hc unused

%Compute M and d
xeq2=xeq(2);
ueq2=ueq(2);
umax2=u_max(2);
umin2=u_min(2);


R2=1;
Q2=1;

ymax2=100; %dummy
ymin2=-100; %dummy
[M2,d2,Qbar2,Rbar2,xeq2_hat,ueq2_hat]=buildMdEq(Cc2,xeq2,ueq2,Q2,R2,ymax2,ymin2,umax2,umin2,ny,nu,pc2);

%Form cost function - constant portion
S2=2*Hx2'*Qbar2*Hx2+Rbar2;

xa2=xa(3);

%Prepare for MPC loop

%Initialize global Lagrange multiplier
%Lambda=[0;0;0];

Lambda=[0.1772;-0.4821;-2.2836]; %hot started from convergence value at first MPC iteration, for problem parameters and gains exactly as given

cr=3; %convergence rate

uprev=zeros(m,1);
xnow=zeros(3,1);
u=uprev;

%Main MPC loop

for k=1:simhor
   T1=0;T2=0; 
   converged=0;
   %Lambda=[0;0;0];
   xa1=xa(1:2);
   xa2=xa(3);
   
   while ~converged 
       
       f1=2*((xa1'*Px1'-xeq1_hat')*Qbar1*Hx1-ueq1_hat'*Rbar1)+Lambda'*R1c;

       f2=2*((xa2'*Px2'-xeq2_hat')*Qbar2*Hx2-ueq2_hat'*Rbar2)+Lambda'*R2c;

       t1=tic;
       Y1=quadprog(S1,f1,M1,d1); %Here we handle inequality constraints at the local level
       %Dual decomposition is used only for the terminal equality
       %constraints 
       
       dt1=toc(t1);
       T1=T1+dt1;
       
       t2=tic;
       Y2=quadprog(S2,f2,M2,d2); 
       dt2=toc(t2);
       T2=T2+dt2;

       %Find constraint error

       %interleave Y1 and Y2 solutions HARD-CODED FOR GIVEN DIMENSIONS AND
       %HORIZON 
       
       Y=[Y1(1);Y2(1);Y1(2);Y2(2);Y1(3);Y2(3);Y1(4);Y2(4);Y1(5);Y2(5)];
       ec=Qc*xa+Aeq*Y-xeq;

       %Central Lagrangian update
       Lambda=Lambda+cr*ec; 
       normEC=norm(ec)
       converged=normEC<0.01;

   end
   
   Ttotal1(k)=T1; %total time spent by Lagrangian iterations, for each "processor"
   Ttotal2(k)=T2;
   
   u_apply=[Y1(1);Y2(1)]; %extract first term of optimal sequence
   u=[u u_apply];  %store control history
   xnext(:,k)=Ad*xnow+Bd*u(:,end); %update plant
   xnow=xnext(:,k); 
   uprev=u(:,end);
   xa=xnow;
end

figure(1)
subplot(3,1,1)
stairs(xnext(1,:),'r--')
subplot(3,1,2)
stairs(xnext(2,:),'r--')
subplot(3,1,3)
stairs(xnext(3,:),'r--')
ylabel('x_3', 'Fontsize',14)
xlabel('MPC iteration', 'Fontsize',14)
legend('Centralized','Distrib')
subplot(3,1,1)
title('States: Centralized vs. Parallel MPC','Fontsize',14)
ylabel('x_1', 'Fontsize',14)
legend('Centralized','Distrib')
subplot(3,1,2)
ylabel('x_2', 'Fontsize',14)
legend('Centralized','Distrib')
subplot(3,1,3)
ylabel('x_3', 'Fontsize',14)
legend('Centralized','Distrib')

figure(2);
subplot(2,1,1)
stairs(u(1,:),'r--')
subplot(2,1,2)
stairs(u(2,:),'r--')
ylabel('u_2', 'Fontsize',14)
xlabel('MPC iteration', 'Fontsize',14)
legend('Centralized','Distrib')
subplot(2,1,1)
title('Control Inputs: Centralized vs. Parallel MPC - Input Constraints','Fontsize',14)
ylabel('u_1', 'Fontsize',14)
legend('Centralized','Distrib')


figure(3)
stairs(Ttotal1,'r--');hold on
stairs(Ttotal2,'k--');
stairs(T,'k')
title('Computation Time: Centralized vs. Parallel MPC - Input Constraints','Fontsize',14)
ylabel('Time per quadratic optimization', 'Fontsize',14)
xlabel('MPC iteration', 'Fontsize',14)
legend('Distrib, processor 1','Distrib, processor 2','Centralized quadprog')

end

function [Px,Hx]=buildPH(A,B,nu,ny,p)

n=size(A,1);m=size(B,2);

%Compute Px
Px=A; 
for i=1:ny-1
  Px=[Px;A^(i+1)];
end

%Compute Hx
Hx=zeros(p*ny,m*nu);
for i=1:ny,
    for j=1:i,
        Hx(1+(i-1)*p:i*p,1+(j-1)*m:j*m)=A^(i-j)*B;
    end
end

end


function [Pc,Hc,Cc]=buildPHC(Cconstr,A,B,ny,nu,pc)

m=size(B,2);

Pc=Cconstr*A;
for i=1:ny-1,
  Pc=[Pc;Cconstr*A^(i+1)];
end

%Compute H for constraints

Hc=zeros(pc*ny,m*nu);
for i=1:ny,
    for j=  1:i,
        Hc(1+(i-1)*pc:i*pc,1+(j-1)*m:j*m)=Cconstr*A^(i-j)*B;
    end
end

%Compute Cc 
Cc=zeros(m*nu,m*nu);
for i=1:nu,
    for j=1:i,
        Cc(1+(i-1)*m:i*m,1+(j-1)*m:j*m)=eye(m);
    end
end

end

function [M,d,Qbar,Rbar,xeq_hat,ueq_hat]=buildMdEq(Cc,xeq,ueq,Q,R,ymax,ymin,umax,umin,ny,nu,pc)


Ly=eye(pc);
xeq_hat=xeq;
ueq_hat=ueq;
Qbar=Q;
Rbar=R;
m=size(R,1);

for i=1:ny-1,
    Ly=[Ly;eye(pc)];
    xeq_hat=[xeq_hat;xeq];
    ueq_hat=[ueq_hat;ueq];
    Qbar=blkdiag(Qbar,Q);
    Rbar=blkdiag(Rbar,R);
end

y_max=Ly*ymax;
y_min=Ly*ymin;

Lu=eye(m);
for i=1:nu-1,
    Lu=[Lu;eye(m)];
end

dup=Lu*umax;
dum=-Lu*umin;

M=[Cc;-Cc];
d=[dup;dum];

end


