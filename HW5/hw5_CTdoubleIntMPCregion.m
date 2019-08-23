
function hw5_CTdoubleIntMPCregion(P)

% clear all
% close all
clc

% Set options 
tol_opt       = 1e-3;
options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 50000,...
                'FinDiffType', 'forward',...
                'RelLineSrchBnd', [],...
                'RelLineSrchBndDuration', 1,...
                'TolConSQP', 1e-6);

n = 2; % state dimension
m = 1; % input dimension 

% eqilibilium point
x_eq = [0;0];
u_eq = 0;

% Number of MPC iterations
mpciterations = 140; 

% Horizon (continuous)
T = 0.5;

% sampling time
delta = 0.05;

% Horizon
N = floor(T/delta);

% initial conditions
t_0 = 0.0;
xplant = [-0.5; -0.7];

%running cost
Q = eye(n);
R = 1;

Umax=1;

% Initial guess for states by simulation
x0=zeros(n*(N+1),1); %preallocate

% Initial guess for input
u0 = Umax*ones(m*N,1)/2;
x0(1:n) = xplant;
for k=1:N
     x0(n*k+1:n*(k+1)) = dynamics(delta,x0(n*(k-1)+1:n*k), u0(k));
end

% Linear constraints
%Inequality constraints: 
A=[];
b=[];
%initial and final state constraint
Aeq=[eye(n),zeros(n,N*(n+m))];
% Aeq = [eye(n), zeros(n,N*(n+m)); zeros(n,N*n), eye(n), zeros(n,N*m)];

beq=xplant;
%input constraints;
lb=[-inf*ones(n*(N+1),1);-Umax*ones(m*N,1)];
ub=[inf*ones(n*(N+1),1);Umax*ones(m*N,1)];
% Set variables for output
t = [];
x = [];
u = [];

% f1 = figure(1); hold on
%plot initial point
fig5 = figure('Name','Question 5','units','normalized', ...
    'outerposition',[0 0 1 1], 'NumberTitle','off');

% Print Header
fprintf('   k  |      u(k)        x(1)        x(2)     Time \n');
fprintf('---------------------------------------------------\n');

% initilization of measured values
tcloop = t_0;
xcloop = xplant;

% simulation
for ii = 1:mpciterations % maximal number of iterations
    
    
    % Set initial guess and initial constraint
    beq=xcloop;
%     beq = [xcloop; x_eq];
    y_init=[x0;u0];
    
    t_Start = tic;
    
    % Solve optimization problem
    %structure: y_OL=[x_OL,u_OL];
    [y_OL, V, exitflag, output]=fmincon(@(y) costfunction( N, y, x_eq, u_eq, Q, R,n,m,delta, P),...
        y_init,A,b,Aeq,beq,lb,ub,...
        @(y) nonlinearconstraints(N, delta, y, x_eq, n,m), options);
    t_Elapsed = toc( t_Start );  
    
    x_OL=y_OL(1:n*(N+1));
    u_OL=y_OL(n*(N+1)+1:end);
    
    %%    
 
    % Store closed loop data
    t = [ t, tcloop ];
    x = [ x, xcloop ];
    u = [ u, u_OL(1:m) ];
    
    % Update closed-loop system (apply first control move to system)
    xcloop = x_OL(n+1:2*n);
    tcloop = tcloop + delta;
        
    % Compute initial guess for next time step
    u0 = [u_OL(m+1:end); zeros(m,1)];
    x0 = [x_OL(n+1:end); dynamics(delta, x_OL(end-n+1:end), u0(end-m+1:end))];
    %%
    % Print numbers
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f  %+6.3f\n', ii, u(end),...
            x(1,end), x(2,end),t_Elapsed);
    
    %plot predicted and closed-loop state trajetories 
    set(0, 'CurrentFigure', 4)
    ax1 = subplot(121); grid on; hold on;    
    plot(x(1,:),x(2,:),'b'), grid on, hold on,
    plot(x_OL(1:n:n*(N+1)),x_OL(n:n:n*(N+1)),'g')
    plot(x(1,:),x(2,:),'ob')
    xlabel('x(1)')
    ylabel('x(2)')
    strTitle = ['Running MPC with Terminal Region ... #',num2str(ii),' /',num2str(mpciterations)];
    title(strTitle);
    axis equal
    drawnow
  
end
ax2 = subplot(122);  grid on; hold on; yyaxis left
stairs(t,u);

end

function xdot = system(~, x, u)
    %State derivatives
    xdot = zeros(2,1);
    xdot(1) = x(2) ;
    xdot(2) =  u;
    
end

function cost = costfunction(N, y, x_eq, u_eq, Q, R,n,m,delta, P)
    
    cost = 0;
    x=y(1:n*(N+1));
    u=y(n*(N+1)+1:end);
    
    % Build the cost by summing up the stage cost and the
    % terminal cost
    for k=1:N
        x_k=x(n*(k-1)+1:n*k);
        u_k=u(m*(k-1)+1:m*k);
        cost = cost + delta*runningcosts(x_k, u_k, x_eq, u_eq, Q, R) + terminalCost(x_k, x_eq, P);
    end
    
end

function cost = runningcosts(x, u, x_eq, u_eq, Q, R)

    cost = (x-x_eq)'*Q*(x-x_eq) + (u-u_eq)'*R*(u-u_eq);
    
end

function alpha = terminalCost(x, x_eq, P)
    alpha = (x-x_eq)'*P*(x-x_eq);
end

  function [c, ceq] = nonlinearconstraints(N, delta, y, x_eq,n,m) 
     
   x=y(1:n*(N+1));
   u=y(n*(N+1)+1:end);
   c = [];
   ceq = [];
 
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);        
        u_k=u((k-1)*m+1:k*m);
        %dynamic constraint
        ceqnew=x_new - dynamics(delta, x_k, u_k);
        ceq = [ceq ceqnew];
    end
  end


function [x] = dynamics(delta, x0, u)
    %integrate numerically
    
    % Set options
    atol_ode  = 1e-4;
    rtol_ode  = 1e-4;
    options = odeset('AbsTol', atol_ode, 'RelTol', rtol_ode);
    
    [t_end, x_end] = ode45(@system, [0,delta], x0, options, u);
        
    x = x_end(end,:)';

end