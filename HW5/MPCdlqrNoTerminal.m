
%Example: MPC for a discrete-time system with LQR cost and no
%terminal conditions
%Multiple-shooting implementation with fmincon

function [t, x, u] = MPCdlqrNoTerminal
% clear;close all;clc

% Set options 
tol      = 1e-3;
options = optimset('Display','on',...
                'TolFun', tol,...
                'MaxIter', 10000,...
                'FinDiffType', 'forward',...
                'RelLineSrchBnd', [],...
                'RelLineSrchBndDuration', 1,...
                'TolConSQP', 1e-6);

% system dimensions
n = 2; %state dimension
m = 1; %input dimension 

%terminal state
x1eq=0;
x2eq=0;

x_eq = [x1eq;x2eq];
u_eq = 0;

%Simulation horizon
simhor = 20;

%Optimization Horizon
N = 3;

% initial conditions for the plant 
t0 = 0.0;
x0plant = [2;2];

% Initial control guess
u0 = zeros(m*N,1);
% Initial state trajectory guess
x0=zeros(n*(N+1),1);

% Linear constraints
%Inequality constraints: 
A=[];
b=[]; %Use these to include state and output constraints if applicable

%initial state equality constraint
Aeq=[eye(n),zeros(n,N*(n+m))];

%beq computed inside mpc loop

%input constraints
umax=20;
lb=[-inf*ones(n*(N+1),1);-umax*ones(m*N,1)];
ub=[+inf*ones(n*(N+1),1);+umax*ones(m*N,1)];

%Cost weights
Q=eye(2);
R=1;

% Initialize data logging

t = [];
x = [];
u = [];

% initilization of closed-loop state and time
tcloop = t0;
xcloop = x0plant;

%plot initial point
% figure(1)
% plot(tcloop,xcloop(1),'bo'), grid on, hold on,
% plot(tcloop,xcloop(2),'r*')
ax5 = subplot(324); grid on; hold on; xlim([0 simhor])
ax4 = subplot(322); grid on; hold on; xlim([0 simhor])

plot(tcloop,xcloop(1),'bo')
plot(tcloop,xcloop(2),'r*')

legend('off')

% simulation
for i = 1:simhor % maximal number of iterations
    
    % Update initial state equality constraint
    beq=xcloop;
      
    % Set initial guess and initial constraint
    y0=[x0;u0];
    
    t1 = tic;
    
    % Solve optimization problem
    
    [y, V, exitflag, output]=fmincon(@(y) cost(N,y,x_eq,u_eq,Q,R,n,m),...
        y0,A,b,Aeq,beq,lb,ub,...
        @(y) nonlinearconstraints(N,y,n,m), options);
    DeltaT = toc( t1 ); 
    
    xpred=y(1:n*(N+1));
    upred=y(n*(N+1)+1:end);
    
    
        
    %Log closed loop trajectories
    t = [ t, tcloop ];
    x = [ x, xcloop ];
    u = [ u, upred(1:m) ];
    
    %Update closed-loop system
    %Note: This example has perfect plant knowledge and no disturbance,
    %so the plant will move as predicted.
    
    xcloop = xpred(n+1:2*n);
    tcloop = tcloop + 1;
        
    %Use previous solution for next initial guess
    u0 = [upred(m+1:end); zeros(m,1)];
    x0 = [xpred(n+1:end); dynamics(xpred(end-n+1:end),u0(end-m+1:end))];
    
    %Show elapsed time
    DeltaT   
    
%     figure(1)
%     plot(tcloop,xcloop(1),'bo'), grid on, hold on,
%     plot(tcloop,xcloop(2),'r*')
%     xlabel('Time')
%     ylabel('x_1,x_2')
%     drawnow
%     figure(2)
%     plot(tcloop-1,u(end),'k+'), grid on, hold on; 
%     drawnow
%     
%     %Plot predicted phase trajectories
%     figure(3)
%     plot(xpred(1:n:n*(N+1)),xpred(n:n:n*(N+1)),'k-*');hold on
    ax4 = subplot(322); grid on; hold on;
    plot(ax4, tcloop,xcloop(1),'bo')
    plot(ax4, tcloop,xcloop(2),'r*')
    xlabel('Time'); ylabel('x_1,x_2');
    strTitle = ['Running MPC with Free Constraints ... #',num2str(i)];
    title(strTitle);
    drawnow
    
    ax5 = subplot(324); grid on; hold on;
    plot(tcloop-1,u(end),'k+')
    ylabel('Control Input');
    drawnow
    
    %Plot predicted phase trajectories
    ax6 = subplot(326); grid on; hold on;
    p = plot(xpred(1:n:n*(N+1)),xpred(n:n:n*(N+1)),'k-', 'LineWidth', 2);hold on
    p.Color(4) = 0.1;    xs = xpred(1:n:n*(N+1));
    ys = xpred(n:n:n*(N+1));
    hold on
    scatter1 = scatter(xs,ys,'MarkerFaceColor','k','MarkerEdgeColor','k'); 
    scatter1.MarkerFaceAlpha = .2;
    scatter1.MarkerFaceAlpha = .2;
    hold off
    xlabel('x1'); ylabel('x2');
    title('Phase Plane Predicted Trajectories');
    axis equal
end
axes(ax4);
title('Closed-loop Trajectories with Free Constraints');
legend('x1', 'x2');
end

function [x] = dynamics(x0, u)    
    %One step ahead prediction of plant states 
    x1=-x0(1)+2*x0(2);
    x2=-x0(2)+2*u;
    x=[x1;x2];
end

function cost = cost(N,y,x_eq,u_eq,Q,R,n,m)
    cost = 0;
    u=y(n*(N+1)+1:end); %extract control section of y  
    x=y(1:n*(N+1));     %extract state section of y  
    for k=1:N
        x_k=x(n*(k-1)+1:n*k); %extract state vector at time k
        u_k=u(m*(k-1)+1:m*k); %extract control vector at time k
        cost = cost + runningcosts(x_k,u_k,x_eq,u_eq,Q,R);
    end
end

function l = runningcosts(x,u,x_eq,u_eq,Q,R)
    l =(x-x_eq)'*Q*(x-x_eq) + (u-u_eq)'*R*(u-u_eq);
end

function [c, ceq] = nonlinearconstraints(N,y,n,m) 
   x=y(1:n*(N+1));
   u=y(n*(N+1)+1:end);
   c = [];
   ceq = [];
   % constraints along prediction horizon
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);        
        u_k=u((k-1)*m+1:k*m);
        
        %dynamic constraint
        ceqnew=x_new - dynamics(x_k, u_k);
        ceq = [ceq ceqnew];
    end
end

