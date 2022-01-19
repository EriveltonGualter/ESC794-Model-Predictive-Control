%Elevator case study
clear all;clc;close all
global ny nu

%% System parameters and constraints 
Jt = 20; %total inertia reflected to linear coordinate, kg-m^2
a  = 12; %gear ratio times torque constant, N-m/A
Rm = 0.1; %motor resistance, Ohm
R  = 1; %drum radius, m
mc = 1000; %nominal elevator car mass, with passengers, kg
mw = 750; %counterweight mass, kg
g  = 9.81; %gravity, m/s^2

mcw=mc-mw;

%Specifications
% to achieve the time limit of 12 seconds I had to change the limitations
% on voltage, acceleration and velocity
voltbar = 55;  %max voltage magnitude
Abar = 0.25*g; %max acceleration magnitude, m/s^2
Vbar = 4; %max velocity magnitude, m/s

%Feasibility of steady voltage
Vss=Rm*R*mcw*g/a; %check if the steady voltage to hold the elevator is within constraints
%Vss = 20.43 it is within constraints 

u_max = (a*R*voltbar/Rm)-R^2*mcw*g; %compute upper constraint on u
u_min = -(a*R*voltbar/Rm)-R^2*mcw*g; %compute lower constraint on u
 
%% Continuous-time state-space matrices with pos - vel states
A=[0 1;0 -a^2/(Jt*Rm)];  
B=[0 1/Jt]';

%C, D matrices for various outputs :
Cpos=[1 0]; Dpos=[0];
Cvel=[0 1]; Dvel=[0];
Caccel=[0 -a^2/(Jt*Rm)]; Daccel=[1/Jt];
    
%Discretize with zoh 

Ts=0.05; %sampling period, s
sysc=ss(A,B,Cpos,Dpos);  % ss(A,B,Cvel,Dvel) ss(A,B,Caccel,Daccel)
sysd=c2d(sysc,Ts);
[Ad,Bd,Cd,Dd]=ssdata(sysd);

%Augment the system for position output 
m = 1; n =2;
Ada=[Ad Bd;zeros(m,n) eye(m)];
Bda=[Bd;eye(m)];

%Rename variables for code re-utilization (similar to Richter's book)
Adu=Ad;
Bdu=Bd;
Ad=Ada;
Bd=Bda;

%% Prediction and control horizons and lambda

lambda=0.000001;
lambdau=0.000001;
%Horizons
ny = 10;
nu = 9;
%%
%Prediction matrices for position (to form cost function) 
C=[1 0 0];D=[0];

n=size(Ad,1); %Augmented states dimension
m=size(Bd,2); %Augmented input dimension
p=size(C,1);  %Augmented output dimension
%Compute Ppos
P=C*Ad; 
for i=1:ny-1
  P=[P;C*Ad^(i+1)];
end
Ppos=P;

%Compute Hpos
H=zeros(p*ny,m*nu);
for i=1:ny
    for j=1:i
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m)=C*Ad^(i-j)*Bd;  
    end
    H(1+(i-1)*p:i*p,1+(j-1)*m+m:j*m+m)=D;
end
%Retain only the first nu blocks (control horizon shorter than prediction horizon)
H=H(:,1:nu*m);
Hpos=H;
%%
%Compute Pvel
C=[0 1 0];D=[0];
P=C*Ad; 
for i=1:ny-1
  P=[P;C*Ad^(i+1)];
end
Pvel=P;

%Compute Hvel
H=zeros(p*ny,m*nu);
for i=1:ny
    for j=1:i
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m)=C*Ad^(i-j)*Bd;  
    end
    H(1+(i-1)*p:i*p,1+(j-1)*m+m:j*m+m)=D;
end
%Retain only the first nu blocks (control horizon shorter than prediction horizon)
H=H(:,1:nu*m);
Hvel=H;

%%
%Prediction matrices for control u
C=[0 0 1];D=[1];  %zeros(m,n),eye(m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p=size(C,1); %output dimension

%Compute Pu
%P=[];
P=C*Ad; 
for i=1:ny-1
  P=[P;C*Ad^(i+1)];
end
Pu=P;
%Compute Hu
%H=[];
H=zeros(p*ny,m*nu);
for i=1:ny
    for j=1:i
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m)=C*Ad^(i-j)*Bd;
    end
    H(1+(i-1)*p:i*p,1+(j-1)*m+m:j*m+m)=D;
end
%Retain only the first nu blocks
H=H(:,1:nu*m);
Hu=H;

%Matrices for constraints on velocity, acceleration
C=[0 1 0;0 -a^2/(Rm*Jt) 1/Jt];D=[0 1/Jt]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pc=size(C,1); %output dimension
p=pc;
%Compute Pc (constraints on velocity and acceleration)
P=C*Ad; 
for i=1:ny-1
  P=[P;C*Ad^(i+1)];
end
Pc=P;
%Compute Hc (constraints on v and a)
H=zeros(p*ny,m*nu);
for i=1:ny
    for j=1:i
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m)=C*Ad^(i-j)*Bd;
    end
    H(1+(i-1)*p:i*p,1+(j-1)*m+m:j*m+m)=D;
end
%Retain only the first nu blocks
H=H(:,1:nu*m);
Hc=H;

%Matrices for voltage constraints (taken as constraints on u)
%Compute Cc 
Cc=zeros(m*nu,m*nu);
for i=1:nu
    for j=1:i
        Cc(1+(i-1)*m:i*m,1+(j-1)*m:j*m)=eye(m);
    end
end

%Weights

pos_target=30; %target elevator position, m 
yf=pos_target*ones(ny,1);

%Form the L vectors for output and input constraints
Ly=repmat(eye(pc),[ny 1]);

Ybar=[Vbar;Abar]; % y here means output which is the velocity and acceleration
y_max=Ly*Ybar;
y_min=-y_max;

Lu=repmat(eye(m),[nu 1]);

%Rate limits
dUMAX=1000; %V/s
UB=dUMAX*ones(nu,1);
LB=-UB;


%D-T Simulation: initial position assumed to be equilibrium at zero
x=[0;0]; %Note that u=0 at the initial time, but V is not zero
xa=zeros(n,1);

%**************************************************************************
RMG = repmat(mcw*g,[ny 1]);

S=Hpos'*Hpos+ (lambda*Rm/(a^2*R^2))*Hu'*Hu+lambdau*eye(m*nu) -(lambda/R^2)*Hu'*Hvel;
S = 2*S;


f=2*(xa'*Ppos'-yf')*Hpos+ 2*(lambda*Rm/(a^2*R^2))*(xa'*Pu')*Hu -(lambda/R^2)*xa'*(Pu'*Hvel+Pvel'*Hu) + (2*lambda*g*Rm/a^2)*RMG'*Hu - lambda*RMG'*Hvel; 
f = f';
%**************************************************************************
dyp=y_max-Pc*xa;
dym=-y_min+Pc*xa;

uprev=zeros(m,1); %start at equilibrium
dup=Lu*(u_max-uprev);
dum=-Lu*(u_min-uprev);

M=[Hc;-Hc;Cc;-Cc];
d=[dyp;dym;dup;dum];

xnow=x; %unaugmented state
ynext=[Cpos;Cvel;Caccel]*xnow+[Dpos;Dvel;Daccel]*uprev;
u=uprev;
du=0;
%Enter MPC loop
simhor=240; %simulation horizon in seconds
for k=1:simhor
   ctrlvec=quadprog(S,f,M,d,[],[],LB,UB);
   u_apply=ctrlvec(1:m); %extract first term of optimal sequence this is del_u
   u=[u uprev+u_apply];  %store control history
   du=[du;u_apply];  %store incremental control history
   xnext(:,k)=Adu*xnow+Bdu*u(:,end); %update plant
   ynext=[ynext [Cpos;Cvel;Caccel]*xnext(:,k)+[Dpos;Dvel;Daccel]*u(:,end)]; %store output history
   xnow=xnext(:,k); 
   uprev=u(:,end);
   xa=[xnow;uprev];

   f=2*(xa'*Ppos'-yf')*Hpos+ 2*(lambda*Rm/(a^2*R^2))*(xa'*Pu')*Hu -(lambda/R^2)*xa'*(Pu'*Hvel+Pvel'*Hu) + (2*lambda*g*Rm/a^2)*RMG'*Hu - lambda*RMG'*Hvel; 
   f = f';
   
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

