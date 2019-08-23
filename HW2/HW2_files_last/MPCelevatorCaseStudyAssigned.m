%Elevator case study
clear;clc;close all

%System parameters
Jt=; %total inertia reflected to linear coordinate, kg-m^2
a=; %gear ratio times torque constant, N-m/A
Rm=; %motor resistance, Ohm
R=; %drum radius, m
mc=; %nominal elevator car mass, with passengers, kg
mw=; %counterweight mass, kg
g=; %gravity, m/s^2

mcw=;

%Specifications
voltbar=;  %max voltage magnitude
Abar=; %max acceleration magnitude, m/s^2
Vbar=; %max velocity magnitude, m/s

%Feasibility of steady voltage
Vss= %check if the steady voltage to hold the elevator is within constraints

u_max=; %compute upper constraint on u
u_min=; %compute lower constraint on u


%Continuous-time state-space matrices with pos - vel states
A=;  B=;

%C, D matrices for various outputs :
Cpos=; Dpos=;
Cvel=; Dvel=;
Caccel=; Daccel=;
    
%Discretize with zoh (check that it preserves state meanings)

Ts=; %sampling period, s
sysc=;
sysd=;
[Ad,Bd,Cd,Dd]=ssdata(sysd);

%Augment the system
Ada=;
Bda=;

%Rename variables for code re-utilization (similar to Richter's book)
Adu=Ad;
Bdu=Bd;
Ad=Ada;
Bd=Bda;


%Prediction and control horizons
%Horizons
ny=;nu=;


%Prediction matrices for position (to form cost function) 
C=;D=;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=size(Ad,1); %state dimension
m=size(Bd,2); %input dimension
p=size(C,1); %output dimension

%Compute Ppos
P=[];
P=; 
for i=1:ny-1
  P=;
end
Ppos=P;
%Compute Hpos
H=[];
H=zeros(p*ny,m*nu);
for i=1:ny,
    for j=1:i,
        H()=C*Ad^(i-j)*Bd;
    end
    H()=D;
end
%Retain only the first nu blocks (control horizon shorter than prediction horizon)
H=H(:,1: );
Hpos=H;

%Prediction matrices for control u
C=[];D=;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=size(C,1); %output dimension

%Compute Pu
P=[];
P=; 
for i=1:ny-1,
  P=[P;];
end
Pu=P;
%Compute Hu
H=[];
H=zeros(p*ny,m*nu);
for i=1:ny,
    for j=1:i,
        H(:)=C*Ad^(i-j)*Bd;
    end
    H(:i)=;
end
%Retain only the first nu blocks
H=H(:,1:);
Hu=H;


%Matrices for constraints on velocity, acceleration
C=[;];D=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pc=size(C,1); %output dimension
p=pc;
%Compute Pc (constraints on v and a)
P=[];
P=C*Ad; 
for i=1:ny-1,
  P=[];
end
Pc=P;
%Compute Hc (constraints on v and a)
H=[];
H=zeros();
for i=1:ny,
    for j=1:i,
        H()=C*Ad^(i-j)*Bd;
    end
    H(:)=D;
end
%Retain only the first nu blocks
H=H(:,1:);
Hc=H;

%Matrices for voltage constraints (taken as constraints on u)
%Compute Cc 
Cc=zeros(m*nu,m*nu);
for i=1:nu,
    for j=1:i,
        Cc(1+(i-1)*m:i*m,1+(j-1)*m:j*m)=;
    end
end


%Weights
lambda=;
lambdau=;

pos_target=; %target elevator position, m 
yf=pos_target*ones(,);

%Form the L vectors for output and input constraints
Ly=repmat(eye(pc),[ny 1]);

Ybar=[Vbar;Abar];
y_max=Ly*Ybar;
y_min=-y_max;

Lu=repmat(,);

%Rate limits
dUMAX=1000; %V/s
UB=dUMAX*ones(nu,1);
LB=-UB;

simhor=250; %simulation horizon in seconds
S=Hpos'*Hpos+Hu'*Hu+lambdau*eye(m*nu); %Hessian portion (constant)

%D-T Simulation: initial position assumed to be equilibrium at zero
x=[0;0]; %Note that u=0 at the initial time, but V is not zero
xa=zeros(n,1);


%Initialize
f=2*(Hpos'*(Ppos*xa-yf)+(lambda*Rm/a^2)*Hu'*(Pu*xa/R^2+g*mcw));
dyp=y_max-Pc*xa;
dym=-y_min+Pc*xa;

uprev=0; %start at equilibrium
dup=Lu*(u_max-uprev);
dum=-Lu*();

M=[];
d=[];

xnow=x; %unaugmented state
ynext=[Cpos;Cvel;Caccel]*xnow+[Dpos;Dvel;Daccel]*uprev;
u=uprev;
du=0;
%Enter MPC loop
for k=1:simhor;
   ctrlvec=quadprog(,,,[],[],LB,UB);
   u_apply=ctrlvec(1:m); %extract first term of optimal sequence
   u=[u uprev+u_apply];  %store control history
   du=[du;u_apply];  %store incremental control history
   xnext(:,k)=Adu*xnow+Bdu*u(:,end); %update plant
   ynext=[ynext [Cpos;Cvel;Caccel]*xnext(:,k)+[Dpos;Dvel;Daccel]*u(:,end)]; %store output history
   xnow=xnext(:,k); 
   uprev=u(:,end);
   xa=[xnow;uprev];
   f=2*(Hpos'*(Ppos*xa-yf)+(lambda*Rm/a^2)*Hu'*(Pu*xa/R^2+g*mcw));%reformulate MPC problem
   dyp=y_max-Pc*xa;
   dym=-y_min+Pc*xa;
   dup=Lu*(u_max-uprev);
   dum=-Lu*(u_min-uprev);
   d=[dyp;dym;dup;dum];
end

%Calculate voltage and current histories
V=Rm*(u+R^2*mcw*g)/a/R;
I=(V-a*ynext(2,:)/R)/Rm;
%Power and energy
Pow=I.*V;
Energy=Ts*trapz(Pow)

%Create time vector for plotting
t=Ts*[0:simhor];
plot(t,ynext(1,:),'k')
hold on
plot(t,ynext(2,:),'r')
title('Outputs')
legend('y_1','y_2')

figure(2)
plot(t,V,'k')

figure(3)
plot(t,Pow)
