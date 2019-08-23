%escSDMPCTraffic.m
%Serial distributed MPC for a traffic intersection with autonomous vehicles

%Programmed by H Richter for MPC course: Selected Topics on Engineering
%Science

%Cleveland State University, Mechanical Engineering Department
%Fall 2018


function [car1, car2, car3, car4, tsim] = escSDMPCTraffic
    % close all
    % clearvars -global

    % Set options 
    tol_opt       = 1e-3;
    options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 50000,...
                'Algorithm', 'active-set',...
                'FinDiffType', 'forward',...
                'RelLineSrchBnd', [],...
                'RelLineSrchBndDuration', 1,...
                'TolConSQP', 1e-6);

    % system dimensions
    n = 8; % state dimension
    m = 4; % input dimension 


    %Cars 1 and 2 on x axis, initially positioned at opposite sides of intersection
    %Cars 3 and 4 on y axis, initially positioned at opposite sides of intersection
    %1 and 3 on positive semiaxes

    %State vector
    %[q1 q2 q3 q4 q1dot q2dot q3dot q4dot]';

    %target equilibrium point (at opposite sides) 
    q1eq=-10;
    q2eq= 10;
    q3eq= -10;
    q4eq= 10;

    u1eq=0;
    u2eq=0;
    u3eq=0;
    u4eq=0;

    x_eq = [q1eq;0;q2eq;0;q3eq;0;q4eq;0];
    u_eq = [u1eq;u2eq;u3eq;u4eq];

    % Number of MPC iterations
    mpciterations = 30;

    %Continuous horizon
    T = 1; 

    % sampling time 
    delta = 0.2;

    % Horizon (discrete)
    N = floor(T/delta);

    % initial conditions
    t_0 = 0.0;
    q10=10;
    q20=-10;
    q30=10;
    q40=-10;
    x_init = [q10;0;q20;0;q30;0;q40;0];

    % stage cost
    Q = 0.5*eye(n);

    % Initial guess for input
    UMAX=1000;


    u0 = zeros(m*N,1); %These are the forces applied to the double integrators

    %Better guess:
    N12=floor(N/2);
    N14=floor(N/4);
    N4=N-N12-N14;

    Ug12=2*abs(q1eq)/N14^2/delta^2; %These are the max accelerations to be applied so that the terminal point is reached
    Ug34=2*abs(q3eq)/N12^2/delta^2; %this is based on the physics formula xf=x0+vo*t+(1/2)acc*t^2 adjusted for discretization

    u01=[zeros(N12,1);-Ug12*ones(N14,1);Ug12*ones(N4,1)]; %head start, then all-out
    u02=-u01;
    u03=[-Ug34*ones(N12,1);Ug34*ones(N-N12,1)];  %start immediately
    u04=-u03;

    u0=[];
    for i=1:N,
    u0=[u0;u01(i);u02(i);u03(i);u04(i)];
    end


    % Initial guess for states by simulation
    x0=zeros(n*(N+1),1); %preallocate
    x0(1:n) = x_init;
    for k=1:N
        x0(n*k+1:n*(k+1)) = dynamics(delta,x0(n*(k-1)+1:n*k), u0(m*(k-1)+1:m*k));
    end

    %Separate for plotting
    x01=[];x02=[];x03=[];x04=[];
    for k=1:N
    x01=[x01;x0(n*(k-1)+1)];
    x02=[x02;x0(n*(k-1)+3)];
    x03=[x03;x0(n*(k-1)+5)];
    x04=[x04;x0(n*(k-1)+7)];
    end

    plot(sqrt(x01.^2+x03.^2),'ko');hold on;
    plot(sqrt(x01.^2+x04.^2),'ro');
    plot(sqrt(x02.^2+x03.^2),'k+');
    plot(sqrt(x02.^2+x04.^2),'r+');


    % Linear constraints
    %Inequality constraints: 
    A=[];
    b=[];
    %initial and terminal state equality constraints
    Aeq=[eye(n),zeros(n,N*(n+m));zeros(n,N*n) eye(n) zeros(n,N*m)];
    %beq computed inside mpc loop
    lb=[];
    ub=[];

    % Set variables for output
    t = [];
    x = [];
    a = [];
    u=[];

    % initilization of measured values
    tmeasure = t_0;
    xmeasure = x_init;
    figure(2)
    plot(xmeasure(1),1,'k*');hold on
    plot(xmeasure(3),-1,'k*');
    plot(1,xmeasure(5),'r*');
    plot(-1,xmeasure(7),'r*');

    car1 = [xmeasure(1); 1];
    car2 = [xmeasure(3); -1];
    car3 = [1; xmeasure(5)];
    car4 = [-1; xmeasure(7)];

    for ii = 1:mpciterations 

    beq=[xmeasure;x_eq];
    y_init=[x0;u0];% Set initial guess and initial constraint

    t_Start = tic;

    [y_OL, V, exitflag, output]=fmincon(@(y) quadcost( N, y, x_eq, Q,n,delta),...
        y_init,A,b,Aeq,beq,lb,ub,...
        @(y) nonlinearconstraints(N, delta, y,n,m), options);

    t_Elapsed = toc( t_Start ); 

    x_OL=y_OL(1:n*(N+1));
    u_OL=y_OL(n*(N+1)+1:end);


    % Store closed loop data
    t = [ t, tmeasure ];
    x = [ x, xmeasure ];
    u = [ u, u_OL(1:m) ];

    % Update closed-loop system (apply first control move to system)
    xmeasure = x_OL(n+1:2*n);
    tmeasure = tmeasure + delta;

    % Compute initial guess for next time step
    u0 = [u_OL(m+1:end); zeros(m,1)];
    x0 = [x_OL(n+1:end); dynamics(delta, x_OL(end-n+1:end), u0(end-m+1:end))];

    car1 = [car1 [xmeasure(1); 1]];
    car2 = [car2 [xmeasure(3); -1]];
    car3 = [car3 [1; xmeasure(5)]];
    car4 = [car4 [-1; xmeasure(7)]];

    end
    tsim = 0:delta:delta*mpciterations; 
    animate(car1, car2, car3, car4, tsim, 1)
end

function [x] = dynamics(delta, x0, u)
    %discretize the system using a Runge-Kutta scheme (ode45)
    
    % Set options
    atol_ode  = 1e-4;
    rtol_ode  = 1e-4;
    options = odeset('AbsTol', atol_ode, 'RelTol', rtol_ode);
    
    % Evaluate the system dynamics (integration)
    [t_intermediate, x_intermediate] = ode23(@system, [0,delta], x0, options, u);
        
    x = x_intermediate(end,:)';

end

function xdot =system(~, x, u)
    % System dynamics for 4 decoupled double integrators
    % Used in optimization
    xdot = zeros(8,1); %pre-allocate
    xdot=[x(2);u(1);x(4);u(2);x(6);u(3);x(8);u(4)];    
end

function cost = quadcost(N, y, x_eq, Q, n,delta)

    % Formulate the cost function to be minimized: quadratic penalty on error relative to equilibrium

    %y is a stacked vector of states and controls:
    cost = 0;
    x=y(1:n*(N+1)); %predicted states x1(0),x2(0)...xn(0), x1(1),...xn(1),...x1(N),...xn(N)

    % Build the cost by summing up the stage cost 
    for k=1:N
        x_k=x(n*(k-1)+1:n*k);
        cost = cost + delta*runningcosts(x_k, x_eq, Q);
    end
    %No terminal cost, since equilibrium equality constraint will be used.
end

function cost = runningcosts(x, x_eq, Q)
    % Provide the running cost   
    cost = (x-x_eq)'*Q*(x-x_eq);
end

function [c, ceq] = nonlinearconstraints(N, delta, y,n,m) 
    % Introduce the nonlinear constraints
    x=y(1:n*(N+1));
    u=y(n*(N+1)+1:end);
    c = [];
    ceq = [];
    UMAX=1000;MIND=4;
    % constraints along prediction horizon
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);        
        u_k=u((k-1)*m+1:k*m);
        %dynamic constraint
        ceqnew=x_new - dynamics(delta, x_k, u_k);
        ceq = [ceq ceqnew];
        %nonlinear constraints on input : |u1,2|<=UMAX
        %Impose inequality constraints

        cc=[u_k;-u_k]-[UMAX*ones(4,1);UMAX*ones(4,1)];

        %Calculate distances between vehicles
        d13=sqrt(x_k(1)^2+x_k(5)^2);
        d14=sqrt(x_k(1)^2+x_k(7)^2);
        d23=sqrt(x_k(3)^2+x_k(5)^2);
        d24=sqrt(x_k(3)^2+x_k(7)^2);
        cdist=-[d13 d14 d23 d24]'+MIND*ones(4,1);

        c=[c;cc;cdist];
    end
end
