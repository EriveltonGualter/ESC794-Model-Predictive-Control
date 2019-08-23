%Elevator case study
clear; clc;
% close all

if ~ exist('guihw')
    %System parameters
    Jt  = 20;   % total inertia reflected to linear coordinate, kg-m^2
    a   = 12;   % gear ratio times torque constant, N-m/A
    Rm  = 0.1;  % motor resistance, Ohm
    R   = 1;    % drum radius, m
    mc  = 1000; % nominal elevator car mass, with passengers, kg
    mw  = 750;  % counterweight mass, kg
    g   = 9.81; % gravity, m/s^2

    mcw  = (mc - mw);

    % Rate limits
    dUMAX = 1000; %V/s

    % Specifications
    voltbar = 40;   % max voltage magnitude
    Abar = 0.2*g;     % max acceleration magnitude, m/s^2
    Vbar = 3;       % max velocity magnitude, m/s

    % Weights
    lambda = 1e-3;
    lambdau = 1e-6;

    Vmax = voltbar;
    u_max = a*R*Vmax/Rm - R^2*mcw*g;   % compute upper constraint on u
    u_min = -a*R*Vmax/Rm - R^2*mcw*g; % compute lower constraint on u

    % Target
    pos_target = 30; % target elevator position, m 
    
    simhor = 250*3; % simulation horizon in seconds
    
    % Prediction and control horizons
    % Horizons
    ny = 60; nu = ny-1;
end

% Feasibility of steady voltage
Vss = R^2*mcw*g*Rm/(a*R); % check if the steady voltage to hold the elevator is within constraints

% Continuous-time state-space matrices with pos - vel states
A = [0 1; 0 -a^2/(Rm*Jt)];  B = [0; 1/Jt];

% C, D matrices for various outputs :
Cpos = [1 0]; Dpos = 0;
Cvel = [0 1]; Dvel = 0;
Caccel = [0 -a^2/(Rm*Jt)]; Daccel = 1/Jt;
    
% Discretize with zoh (check that it preserves state meanings)
Ts = 0.05; %sampling period, s
sysc = ss(A, B, eye(2), 0);
sysd = c2d(sysc, Ts, 'zoh');
[Ad,Bd,Cd,Dd] = ssdata(sysd);

% Augment the system
% Since size([Ad Bd]) = 2x3:
Ada = [Ad Bd; zeros(1,2) eye(1)]; 
Bda = [Bd; eye(1)];

% Rename variables for code re-utilization (similar to Richter's book)
Adu = Ad;
Bdu = Bd;
Ad = Ada;
Bd = Bda;

% Prediction matrices for position (to form cost function) 
C = [Cpos 0]; D = Dpos;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(Ad,1); % state dimension
m = size(Bd,2); % input dimension
p = size(C,1); % output dimension

% Compute Ppos
P = [];
P = C*Ad; 
for i=1:ny-1
  P = [P; C*Ad^(i+1)];
end
Ppos = P;

% Compute Hpos
H = [];
H = zeros(p*ny,m*nu);
for i=1:ny,
    for j=1:i,
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m) = C*Ad^(i-j)*Bd;
    end
    H(1+(i-1)*p:i*p,1+j*m:(j+1)*m) = D;
end

% Retain only the first nu blocks (control horizon shorter than prediction horizon)
H = H(:,1:nu*m);
Hpos = H;

C = [Cvel 0]; D = Dvel;
% Compute Pvel
P = [];
P = C*Ad; 
for i=1:ny-1
  P = [P; C*Ad^(i+1)];
end
Pvel = P;

% Compute Hvel
H = [];
H = zeros(p*ny,m*nu);
for i=1:ny,
    for j=1:i,
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m) = C*Ad^(i-j)*Bd;
    end
    H(1+(i-1)*p:i*p,1+j*m:(j+1)*m) = D;
end

% Retain only the first nu blocks (control horizon shorter than prediction horizon)
H = H(:,1:nu*m);
Hvel = H;

% Prediction matrices for control u
C = [0 0 1]; D = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = size(C,1); % output dimension

% Compute Pu
P = [];
P = C*Ad; 
for i=1:ny-1,
  P = [P; C*Ad^(i+1)];
end
Pu = P;

% Compute Hu
H = [];
H = zeros(p*ny,m*nu);
for i=1:ny,
    for j=1:i,
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m) = C*Ad^(i-j)*Bd;
    end
    H(1+(i-1)*p:i*p,1+j*m:(j+1)*m) = D;
end
% Retain only the first nu blocks
H = H(:,1:nu*m);
Hu = H;

% Matrices for constraints on velocity, acceleration
C = [Cvel 0; Caccel Daccel]; D = [Dvel; Daccel];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pc = size(C,1); %output dimension
p = pc;

%Compute Pc (constraints on v and a)
P = [];
P = C*Ad; 
for i=1:ny-1,
  P = [P; C*Ad^(i+1)];
end
Pc = P;

% Compute Hc (constraints on v and a)
H = [];
H = zeros(p*ny,m*nu);
for i=1:ny,
    for j=1:i,
        H(1+(i-1)*p:i*p,1+(j-1)*m:j*m) = C*Ad^(i-j)*Bd;
    end
    H(1+(i-1)*p:i*p,1+j*m:(j+1)*m) = D;
end
% Retain only the first nu blocks
H = H(:,1:nu*m);
Hc = H;

% Matrices for voltage constraints (taken as constraints on u)
% Compute Cc 
Cc=zeros(m*nu,m*nu);
for i=1:nu,
    for j=1:i,
        Cc(1+(i-1)*m:i*m,1+(j-1)*m:j*m) = eye(m);
    end
end

yf = pos_target*ones(ny,1);

% Form the L vectors for output and input constraints
Ly = repmat(eye(pc),[ny 1]);
Lu = repmat(eye(m),[nu 1]);

Ybar = [Vbar; Abar];
y_max = Ly*Ybar;
y_min = -y_max;

% Rate limits
UB = dUMAX*ones(nu,1);
LB = -UB;

% D-T Simulation: initial position assumed to be equilibrium at zero
x  = [0;0]; %Note that u=0 at the initial time, but V is not zero
xa = zeros(n,1);

% Initialize
% f=2*(Hpos'*(Ppos*xa-yf)+(lambda*Rm/a^2)*Hu'*(Pu*xa/R^2+g*mcw));
% S = Hpos'*Hpos+Hu'*Hu+lambdau*eye(m*nu); %Hessian portion (constant)

f = 2*(xa'*Ppos'-yf')*Hpos + 2*(lambda*Rm/(a*R)^2)*xa'*Pu'*Hu + ...
    2*lambda*mcw'*ones(1,ny)*g*Rm*Hu/a^2 - lambda*xa'*(Pu'*Hvel + Pvel'*Hu)/R^2 - ...
    lambda*mcw'*ones(1,ny)*g*Hvel;
S = 2*(Hpos'*Hpos + (lambda*Rm/(a*R)^2)*Hu'*Hu - lambda*Hu'*Hvel/R^2 + lambdau*eye(m*nu));

dyp =  y_max-Pc*xa;
dym = -y_min+Pc*xa;

uprev = 0; %start at equilibrium
dup =  Lu*(u_max-uprev);
dum = -Lu*(u_min-uprev);

M = [Hc; -Hc; Cc; -Cc];
d = [dyp; dym; dup; dum];

xnow = x; %unaugmented state
ynext = [Cpos;Cvel;Caccel]*xnow+[Dpos;Dvel;Daccel]*uprev;
u = uprev;
du = 0;

% Enter MPC loop
for k=1:simhor
   ctrlvec = quadprog(S,f,M,d,[],[],LB,UB);
   u_apply = ctrlvec(1:m); %extract first term of optimal sequence
   u = [u uprev+u_apply];  %store control history
   du = [du;u_apply];  %store incremental control history
   xnext(:,k) = Adu*xnow+Bdu*u(:,end); %update plant
   ynext = [ynext [Cpos;Cvel;Caccel]*xnext(:,k)+[Dpos;Dvel;Daccel]*u(:,end)]; %store output history
   xnow = xnext(:,k); 
   uprev = u(:,end);
   xa = [xnow;uprev];
   
   f=2*(Hpos'*(Ppos*xa-yf)+(lambda*Rm/a^2)*Hu'*(Pu*xa/R^2+g*mcw));
   dyp = y_max-Pc*xa;
   dym = -y_min+Pc*xa;
   dup = Lu*(u_max-uprev);
   dum = -Lu*(u_min-uprev);
   d = [dyp; dym; dup; dum];
end

% Calculate voltage and current histories
V = Rm*(u+R^2*mcw*g)/a/R;
I = (V-a*ynext(2,:)/R)/Rm;

% Power and energy
Pow = I.*V;
Energy = Ts*trapz(Pow)

%Create time vector for plotting
t = Ts*[0:simhor];

ax1 = subplot(221); hold on
plot(t,ynext(1,:),'k')
plot(t,ynext(2,:),'r')
    title('Outputs')
    legend('y_1','y_2')

ax2 = subplot(222); hold on
plot(t,V,'k')
    title('Voltage');

ax3 = subplot(223); hold on
plot(t,Pow)
    title('Power');

ax4 = subplot(224); hold on
plot(t,u)
    title('Control');
