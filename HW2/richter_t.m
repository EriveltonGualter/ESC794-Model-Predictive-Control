%MPC example
%CMAPSS-40k linearized at Ground Idle
%3 constrained inputs
%constraints on 2 outputs: T48 and SM-HPC
%Plant matrices for Ground Idle assumed available
%in the workspace as Aa,Ba,Ca,Da
%Plant matrices

A=Aa;B=Ba;
Cconstr=Ca; D=Da;%to be used for constraints
sysCT=ss(A,B,Cconstr,D);
Ts=0.015;
%Discretize
sysDT=c2d(sysCT,Ts,'zoh');
[Adu,Bdu,Cdu_constr,Ddu]=ssdata(sysDT);
%Augment
Ad=[Adu Bdu;zeros(3,2) eye(3)];Bd=[Bdu;eye(3)];
Cd_constr=[Cdu_constr Ddu];Dd_constr=Ddu;
Cd=[1 0 0 0 0]; %to be used for control
n=size(Ad,1); %state dimension
m=size(Bd,2); %input dimension
p=size(Cd,1); %controlled output dimension
%Horizons
ny=7;nu=3;
%Compute P for control
P=Cd*Ad;
for i=1:ny-1,
P=[P;Cd*Ad^(i+1)];
end
%Compute P for constraints
Pc=Cd_constr*Ad;
for i=1:ny-1,
Pc=[Pc;Cd_constr*Ad^(i+1)];
end
%Compute H for control
H=zeros(p*ny,m*nu);
for i=1:ny,
for j=1:i,
    H(1+(i-1)*p:i*p,1+(j-1)*m:j*m)=Cd*Ad^(i-j)*Bd;
end
end
%Retain only the first nu blocks
H=H(:,1:nu*m);
%Compute H for constraints
pc=size(Cd_constr,1); %constrained output dimension
Hc=zeros(pc*ny,m*nu);
for i=1:ny,
for j=1:i,
Hc(1+(i-1)*pc:i*pc,1+(j-1)*m:j*m)=Cd_constr*Ad^(i-j)*Bd;
end
Hc(1+(i-1)*pc:i*pc,j*m+1:(j+1)*m)=Dd_constr;
end
%Retain only the first nu blocks
Hc=Hc(:,1:nu*m);
%Compute Cc
Cc=zeros(m*nu,m*nu);
for i=1:nu,
for j=1:i,
Cc(1+(i-1)*m:i*m,1+(j-1)*m:j*m)=eye(m);
end
end
%Optimization Weight
lambda=0.01;
r0=100; %controlled output reference
r=100*ones(ny,1);
%Offline computations
Ly=eye(pc);
for i=1:ny-1,
Ly=[Ly;eye(pc)];
end
y_max=Ly*[300;20]; %Incremental upper-bounds: T48 and SmHPC
y_min=Ly*[-150;-10]; %Incremental lower-bounds: T48 and SmHPC
Lu=eye(m);
for i=1:nu-1,
Lu=[Lu;eye(m)];
end
u_max=[2;15;0.4]; %Incremental upper-bounds: WF,VSV and VBV
u_min=[-1;-20;-0.5];%Incremental lower-bounds: WF,VSV and VBV
%Initialization
S=H'*H+lambda*eye(m*nu);
%Initial conditions
xa=zeros(n,1);
f=H'*(P*xa-r);
dyp=y_max-Pc*xa;
dym=-y_min+Pc*xa;
uprev=zeros(m,1);
dup=Lu*(u_max-uprev);
dum=-Lu*(u_min-uprev);
M=[Hc;-Hc;Cc;-Cc];
d=[dyp;dym;dup;dum];
xnow=zeros(2,1); %Plant state initialization
ynext=[1 0;Ca]*xnow;
u=uprev; %Control initialization
simhor=30; %simulation horizon
for k=1:simhor;
ctrlvec=quadprog(S,f,M,d);
u_apply=ctrlvec(1:m); %extract first term of optimal sequence
u=[u uprev+u_apply]; %store control history
xnext(:,k)=Adu*xnow+Bdu*u(:,end); %update plant
ynext=[ynext [1 0;Ca]*xnext(:,k)]; %store output history
xnow=xnext(:,k);
uprev=u(:,end);
xa=[xnow;uprev];
%online MPC reformulation
f=H'*(P*xa-r);
dyp=y_max-Pc*xa;
dym=-y_min+Pc*xa;
dup=Lu*(u_max-uprev);
dum=-Lu*(u_min-uprev);
d=[dyp;dym;dup;dum];
end
t=Ts*[0:simhor];
%Plotting
plot(t,ynext(1,:),'k')
hold on
plot(t,ynext(2,:),'k--')
plot(t,ynext(3,:),'k:')
ylabel('\Delta N_f,\Delta T_{48} and \Delta SmHPC','FontSize',14)
xlabel('Time, sec.', 'FontSize',14)
title('Linearized MPC Simulation: CMPASS-40k Near Ground Idle: Outputs', 'FontSize',14)
legend('\Delta N_f','\Delta T_{48}','\Delta SmHPC')
figure(2)
plot(t,u(1,:),'k')
hold on
plot(t,u(2,:),'k--')
plot(t,u(3,:),'k:')
ylabel('\Delta W_F, \Delta VSV and \Delta VBV','FontSize',14)
xlabel('Time, sec.', 'FontSize',14)
title('Linearized MPC Simulation: CMPASS-40k Near Ground Idle: Inputs', 'FontSize',14)
legend('\Delta W_F','\Delta VSV','\Delta VBV')