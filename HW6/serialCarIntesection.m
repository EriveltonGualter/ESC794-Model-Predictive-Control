% Serial distributed MPC for a traffic intersection with autonomous vehicles

% Erivelton Gualter 
%
% Based on escSDMPCTraffic.m by Dr. Hanz Richter

function [car1, car2, car3, car4, car1s, car2s, car3s, car4s, tsim, T, dt1, dt2, dt3, dt4] = serialCarIntesection
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
%     delta = 0.2;
    delta = 0.1;

    % Horizon (discrete)
    N = floor(T/delta);

    % initial conditions
    t_0 = 0.0;
    q10 = 10;
    q20 = -10;
    q30 = 10;
    q40 = -10;
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

    % Separate for plotting
    x01=[];x02=[];x03=[];x04=[];
    for k=1:N
        x01=[x01;x0(n*(k-1)+1)];
        x02=[x02;x0(n*(k-1)+3)];
        x03=[x03;x0(n*(k-1)+5)];
        x04=[x04;x0(n*(k-1)+7)];
    end

    % Linear constraints
    % Inequality constraints: 
    A=[];
    b=[];
    % Initial and terminal state equality constraints
    Aeq=[eye(n),zeros(n,N*(n+m));zeros(n,N*n) eye(n) zeros(n,N*m)];
    %beq computed inside mpc loop    
    lb = [];
    ub = [];
    
    % Set variables for output
    t = [];
    x = [];

    a = [];
    u = [];
    time_central = 0;
    c = 0;

    % initilization of measured values
    tmeasure = t_0;
    xmeasure = x_init;

    car1 = [xmeasure(1); 1];
    car2 = [xmeasure(3); -1];
    car3 = [1; xmeasure(5)];
    car4 = [-1; xmeasure(7)];

    disp('Finding solution for centralized ...');
    %% Centralized MPC
    for ii = 1:mpciterations 

        beq=[xmeasure;x_eq];
        y_init=[x0;u0];% Set initial guess and initial constraint

        t = tic;
 
        [y_OL, V, exitflag, output]=fmincon(@(y) quadcost( N, y, x_eq, Q,n,delta),...
            y_init,A,b,Aeq,beq,lb,ub,...
            @(y) nonlinearconstraints(N, delta, y,n,m), options);

        T(ii) = toc(t);

        x_OL=y_OL(1:n*(N+1));
        u_OL=y_OL(n*(N+1)+1:end);


        % Store closed loop dataanimate(car1, car2, car3, car4, tsim)

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
        
        % Get initial Condition
        if ii == 1
            % Extracted Initial Condition
            [x01,x02,x03,x04,u01,u02,u03,u04] = getSplitData(y_OL,N,n,m,2,1);
        end
    end
    tsim = 0:delta:delta*mpciterations; 
%     animate(car1, car2, car3, car4, tsim, 1)
    
    %% Serial MPC
    
    time_serial1 = 0;
    time_serial2 = 0;
    time_serial3 = 0;
    time_serial4 = 0;
    
    x_init = [q10;0;q20;0;q30;0;q40;0];

    car1s = [x_init(1); 1];
    car2s = [x_init(3); -1];
    car3s = [1;  x_init(5)];
    car4s = [-1; x_init(7)];
    
    % Initial Condition for subsystems
    x1 = x01;
    u1 = u01;
    x2 = x02;
    u2 = u02;
    x3 = x03;
    u3 = u03;
    x4 = x04;
    u4 = u04;
    
    % initilization of measured values
    tmeasure = t_0;
    xmeasure = x_init;
    
    % Define subsystem dimensions
    ns = 2;
    ms = 1;
    
    % Stage Costs for each subsystem
    Q1 = 1*eye(ns);
    
    % Initial and Terminal state equality constraints
    Aeq = [eye(ns),zeros(ns,N*(ns+ms));zeros(ns,N*ns) eye(ns) zeros(ns,N*ms)];
    
    % Constraint 
    lb=[];
    ub=[];
    
    flags = [];

    disp('Finding solution for serial ...');
    for ii = 1:mpciterations 
        
        if ii == 1
            yneigh = y_OL;
        end
             
        for seq = [1 3 2 4]
            switch seq
                case 1
                    % Substem 1 -------------------------------------------------------
                    Id = 1; carId = 2*Id-1;
                    beq = [xmeasure(carId:carId+1); x_eq(carId:carId+1)];
                    y_init = [x1; u1];       % Set initial guess and initial constraint

                    t1 = tic;
                    
                    [y_OL, V, exitflag, output]=fmincon(@(y) quadcost( N, y, x_eq(carId:carId+1), Q1,ns,delta),...
                        y_init,A,b,Aeq,beq,lb,ub,...
                        @(y) nonlinearconstraints_distrib(N, delta, y, yneigh, n, m, Id) , options);
        
                    dt1(ii) = toc(t1);
                    
                    x01 = y_OL(1:ns*(N+1));
                    u01 = y_OL(ns*(N+1)+1:end);

                    % Initial guess for subsystem 1 optimization: 
                    x1 = [x01(ns+1:end); Decoupleddynamic_distrib(delta, x01(end-ns+1:end), u01(end-ms+1:end))];
                    u1 = [u01(ms+1:end); zeros(ms,1)];
        
                    % Perform communication
                    yneigh = insertNeigh(yneigh, y_OL, N, n, m, ns, ms, Id);
         
                case 2
                    % Substem 2 -------------------------------------------------------
                    Id = 2; carId = 2*Id-1;
                    beq = [xmeasure(carId:carId+1); x_eq(carId:carId+1)];
                    y_init = [x2; u2];       % Set initial guess and initial constraint

                    t2 = tic;
                    
                    [y_OL, V, exitflag, output]=fmincon(@(y) quadcost( N, y, x_eq(carId:carId+1), Q1,ns,delta),...
                        y_init,A,b,Aeq,beq,lb,ub,...
                        @(y) nonlinearconstraints_distrib(N, delta, y, yneigh, n, m, Id) , options);

                    dt2(ii) = toc(t2);
                    
                    x02 = y_OL(1:ns*(N+1));
                    u02 = y_OL(ns*(N+1)+1:end);

                    % Initial guess for subsystem 1 optimization: 
                    x2 = [x02(ns+1:end); Decoupleddynamic_distrib(delta, x02(end-ns+1:end), u02(end-ms+1:end))];
                    u2 = [u02(ms+1:end); zeros(ms,1)];

                    % Perform communication
                    yneigh = insertNeigh(yneigh, y_OL, N, n, m, ns, ms, Id);

                case 3
                    % Substem 3 -------------------------------------------------------
                    Id = 3; carId = 2*Id-1;
                    beq = [xmeasure(carId:carId+1); x_eq(carId:carId+1)];
                    y_init = [x3; u3];       % Set initial guess and initial constraint

                    t3 = tic;
                    
                    [y_OL, V, exitflag, output]=fmincon(@(y) quadcost( N, y, x_eq(carId:carId+1), Q1,ns,delta),...
                        y_init,A,b,Aeq,beq,lb,ub,...
                        @(y) nonlinearconstraints_distrib(N, delta, y, yneigh, n, m, Id) , options);
                    
                    dt3(ii) = toc(t3);
                    
                    x03 = y_OL(1:ns*(N+1));
                    u03 = y_OL(ns*(N+1)+1:end);

                    % Initial guess for subsystem 1 optimization: 
                    x3 = [x03(ns+1:end); Decoupleddynamic_distrib(delta, x03(end-ns+1:end), u03(end-ms+1:end))];
                    u3 = [u03(ms+1:end); zeros(ms,1)];

                    % Perform communication
                    yneigh = insertNeigh(yneigh, y_OL, N, n, m, ns, ms, Id);

                case 4
                    % Substem 4 -------------------------------------------------------
                    Id = 4; carId = 2*Id-1;
                    beq = [xmeasure(carId:carId+1); x_eq(carId:carId+1)];
                    y_init = [x4; u4];       % Set initial guess and initial constraint

                    t4 = tic;
                    
                    [y_OL, V, exitflag, output]=fmincon(@(y) quadcost( N, y, x_eq(carId:carId+1), Q1,ns,delta),...
                        y_init,A,b,Aeq,beq,lb,ub,...
                        @(y) nonlinearconstraints_distrib(N, delta, y, yneigh, n, m, Id) , options);
                    
                    dt4(ii) = toc(t4);

                    x04 = y_OL(1:ns*(N+1));
                    u04 = y_OL(ns*(N+1)+1:end);

                    % Initial guess for subsystem 1 optimization: 
                    x4 = [x04(ns+1:end); Decoupleddynamic_distrib(delta, x04(end-ns+1:end), u04(end-ms+1:end))];
                    u4 = [u04(ms+1:end); zeros(ms,1)];

                    % Perform communication
                    yneigh = insertNeigh(yneigh, y_OL, N, n, m, ns, ms, Id);          
            end
            flags = [flags exitflag]; 
        end
             
        % Store closed loop data
        t = [ t, tmeasure ];
        x = [ x, xmeasure ];

        % Update centralized system with the controls
        u = yneigh((N+1)*n+1:(N+1)*n+m);
        xmeasure = dynamics(delta, xmeasure, u);
        
        tmeasure = tmeasure + delta;   

        car1s = [car1s [xmeasure(1); 1]];
        car2s = [car2s [xmeasure(3); -1]];
        car3s = [car3s [1; xmeasure(5)]];
        car4s = [car4s [-1; xmeasure(7)]];
        
    end
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

%% Split data
function [x01, x02, x03, x04, u1, u2, u3, u4] = getSplitData(X, N, n, m, ns, ms)

    x01(1:ns:(N+1)*ns,1) = X(1:n:(N+1)*n);
    x01(2:ns:(N+1)*ns,1) = X(2:n:(N+1)*n);
    x02(1:ns:(N+1)*ns,1) = X(3:n:(N+1)*n);
    x02(2:ns:(N+1)*ns,1) = X(4:n:(N+1)*n);
    x03(1:ns:(N+1)*ns,1) = X(5:n:(N+1)*n);
    x03(2:ns:(N+1)*ns,1) = X(6:n:(N+1)*n);
    x04(1:ns:(N+1)*ns,1) = X(7:n:(N+1)*n);
    x04(2:ns:(N+1)*ns,1) = X(8:n:(N+1)*n);

    u1 = X((N+1)*n+1:m: end);
    u2 = X((N+1)*n+2:m: end);
    u3 = X((N+1)*n+3:m: end);
    u4 = X((N+1)*n+4:m: end);
end

%% Insert Data
function out = insertNeigh(yneigh, y_O, N, n, m, ns, ms, id)
    
    [x01, x02, x03, x04, u1, u2, u3, u4] = getSplitData(yneigh, N, n, m, ns, ms);

    x_OL = y_O(1:ns*(N+1));
    a_OL = y_O(ns*(N+1)+1:end);

    x01 = (id == 1)*x_OL + (~(id == 1))*x01;
    x02 = (id == 2)*x_OL + (~(id == 2))*x02;
    x03 = (id == 3)*x_OL + (~(id == 3))*x03;
    x04 = (id == 4)*x_OL + (~(id == 4))*x04;
    
    u1 = (id == 1)*a_OL + (~(id == 1))*u1;
    u2 = (id == 2)*a_OL + (~(id == 2))*u2;
    u3 = (id == 3)*a_OL + (~(id == 3))*u3;
    u4 = (id == 4)*a_OL + (~(id == 4))*u4;
        
    out(1:n:(N+1)*n,1) = x01(1:ns:(N+1)*ns);
    out(2:n:(N+1)*n,1) = x01(2:ns:(N+1)*ns);
    out(3:n:(N+1)*n,1) = x02(1:ns:(N+1)*ns);
    out(4:n:(N+1)*n,1) = x02(2:ns:(N+1)*ns);
    out(5:n:(N+1)*n,1) = x03(1:ns:(N+1)*ns);
    out(6:n:(N+1)*n,1) = x03(2:ns:(N+1)*ns);
    out(7:n:(N+1)*n,1) = x04(1:ns:(N+1)*ns);
    out(8:n:(N+1)*n,1) = x04(2:ns:(N+1)*ns);
    
    out((N+1)*n+1:m:(N+1)*n+N*m) = u1;
    out((N+1)*n+2:m:(N+1)*n+N*m) = u2;
    out((N+1)*n+3:m:(N+1)*n+N*m) = u3;
    out((N+1)*n+4:m:(N+1)*n+N*m) = u4;
    
end

%% Nonlinear Constraints
function [c, ceq] = nonlinearconstraints(N, delta, y,n,m) 
    % Introduce the nonlinear constraints
    x=y(1:n*(N+1));
    u=y(n*(N+1)+1:end);
    c = [];
    ceq = [];
    UMAX=1000;MIND=5;
    % constraints along prediction horizon
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);        
        u_k=u((k-1)*m+1:k*m);
        
        % dynamic constraint
        ceqnew = x_new - dynamics(delta, x_k, u_k);
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

% Edit
function [c, ceq] = nonlinearconstraints_distrib(N, delta, yp, yneigh, n, m, id) 
   % Evaluates the nonlinear constraints for a subsystem using neighbor information
   % Introduce the nonlinear constraints
   
   c = [];
   ceq = [];
   
   ns = 2;
   ms = 1;
 
   % Nonlinear Constraints applied on Control Input Neighbor
   UMAX = 1000; MIND = 4;
   
   [x01,x02,x03,x04,u1,u2,u3,u4] = getSplitData(yneigh,N,n,m,ns,ms);
   
   % Dynamic Constraints
   x = yp(1:ns*(N+1));
   u = yp(ns*(N+1)+1:end);
   
   % constraints along prediction horizon
    for k=1:N
        x_k = x((k-1)*ns+1:k*ns);
        x_new = x(k*ns+1:(k+1)*ns);        
        u_k = u((k-1)*ms+1:k*ms);
        
        % Dynamic constraint
        ceqnew = x_new - Decoupleddynamic_distrib(delta, x_k, u_k);
        ceq = [ceq ceqnew];

        xk1 = (id == 1)*x_k(1) + (~(id == 1))*x01((k-1)*ns+1);
        xk2 = (id == 2)*x_k(1) + (~(id == 2))*x02((k-1)*ns+1);
        xk3 = (id == 3)*x_k(1) + (~(id == 3))*x03((k-1)*ns+1);
        xk4 = (id == 4)*x_k(1) + (~(id == 4))*x04((k-1)*ns+1);

        da = 0;
        db = 0;
        if (id == 1) || (id == 2)
            da = x_k(1)^2 + xk3^2;
            db = x_k(1)^2 + xk4^2;
        else
            da = x_k(1)^2 + xk1^2;
            db = x_k(1)^2 + xk2^2;
        end
        
        cc = [-da; -db] + MIND^2;
        
        c = [c; cc];
    end
end

%% Dynamics
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

function [x] = Decoupleddynamic_distrib(delta, x0, a)
    %discretize the system using a Runge-Kutta scheme (ode45)
    
    % Set options
    atol_ode  = 1e-4;
    rtol_ode  = 1e-4;
    options = odeset('AbsTol', atol_ode, 'RelTol', rtol_ode);
    
    % Evaluate the system dynamics (integration)
    [t_intermediate, x_intermediate] = ode23(@Decoupledsystem_distrib, [0,delta], x0, options, a);
        
    x = x_intermediate(end,:)';

end

%% System
function xdot =system(~, x, u)
    % System dynamics for 4 decoupled double integrators
    % Used in optimization
    xdot = zeros(8,1); %pre-allocate
    xdot = [x(2);u(1);x(4);u(2);x(6);u(3);x(8);u(4)];    
end

function xdot = Decoupledsystem_distrib(~, x, a)
    % System dynamics for a double integrator
    xdot=[x(2);a];    
end