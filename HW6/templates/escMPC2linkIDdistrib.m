%MPC2linkID.m
%Inverse dynamics centralized version: the cost function does not include
%control penalties, but joint torque constraints are introduced.

%Input transformation and decoupling through inverse dynamics, optimization
%relative to virtual controls

%Programmed by H Richter for MPC course: Selected Topics on Engineering
%Science

%Numerical solution based on multiple-shooting method - Allgoewer + Mueller

%Cleveland State University, Mechanical Engineering Department
%Fall 2018

function escMPC2linkIDdistrib
    clear all
    close all
    clc

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

    %see what happens when the algorithm is changed from Active Set to default.

    % system dimensions
    n = 4; % state dimension
    m = 2; % input dimension 

    % equilibrium point (computed separately so tau_eq=g(q_eq))
    q1eq=pi/4;
    q2eq=pi/4;
    tau1eq=6.92964645;
    tau2eq=0;
    x_eq = [q1eq;q2eq;0;0];
    u_eq = [tau1eq;tau2eq];

    % Number of MPC iterations
    mpciterations = 10;

    % Horizon (continuous)
    T = 0.5; 

    % sampling time (Discretization steps)
    delta = 0.1;

    % Horizon (discrete)
    N = floor(T/delta);

    % initial conditions
    t_0 = 0.0;
    x_init = [0;0;0;0];

    % stage cost
    Q = 0.5*eye(n);

    % Initial guess for input
    a0 = zeros(m*N,1); %These are virtual accelerations applied to the double integrators
    % Initial gues for states
    x0=zeros(n*(N+1),1);

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
    tau=[];

    % initilization of measured values
    tmeasure = t_0;
    xmeasure = x_init;

    %Initialize by running centralized solution, saving first predictions
    %to be used as guesses for the distributed solution
    subplot(2,2,1)
    plot(t_0,x_init(1),'bo');hold on
    subplot(2,2,3)
    plot(t_0,x_init(2),'ro');hold on

    %% Centralized 
    for ii = 1:mpciterations % maximal number of iterations

        beq=[xmeasure;x_eq];
        y_init=[x0;a0];% Set initial guess and initial constraint

        t_Start = tic;

        % Solve optimization problem
        %structure: y_OL=[x_OL,u_OL];
        [y_OL, V, exitflag, output]=fmincon(@(y) quadcost( N, y, x_eq, Q,n,delta),...
        y_init,A,b,Aeq,beq,lb,ub,...
        @(y) nonlinearconstraints(N, delta, y,n,m), options);

        t_Elapsed = toc( t_Start ); 

        x_OL=y_OL(1:n*(N+1));
        a_OL=y_OL(n*(N+1)+1:end);

        %Calculate actual torque inputs
        taumeasure=computeTau(xmeasure,a_OL(1:m)); 

        % Store closed loop data
        t = [ t, tmeasure ];
        x = [ x, xmeasure ];
        a = [ a, a_OL(1:m) ];
        tau=[tau,taumeasure];
        % Update closed-loop system (apply first control move to system)
        xmeasure = x_OL(n+1:2*n);
        tmeasure = tmeasure + delta;

        % Compute initial guess for next time step
        a0 = [a_OL(m+1:end); zeros(m,1)];
        x0 = [x_OL(n+1:end); Decoupleddynamic(delta, x_OL(end-n+1:end), a0(end-m+1:end))];

        fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f  %+6.3f\n', ii, tau(1,end),...
            tau(2,end),x(1,end),x(2,end),t_Elapsed);
        subplot(2,2,1)
        plot(tmeasure,xmeasure(1),'bo-')
        subplot(2,2,2)
        plot(tmeasure,taumeasure(1),'b*')
        subplot(2,2,3)
        plot(tmeasure,xmeasure(2),'ro-')
        subplot(2,2,4)
        plot(tmeasure,taumeasure(2),'r*')
        xlabel('Time')
        ylabel('Joint Coords')
        drawnow


        if ii==1,       
            [yp,yneigh]=splitdata(y_OL,1,N); %sequence will be 1,2,1,2,...

            %Prepare initial guesses
            % Initial guess for subsystem 1 optimization
            a01 = [a_OL(1);zeros(N-1,1)]; 

            x_eq1=[pi/4;0];
            a_OL2=[yneigh(2*N+4:end);0];
            x_OL2=[yneigh(3:2*(N+1));x_eq1];

            % Initial guess for states: from split data
            x01=[yp(3:2*(N+1));x_eq1];
        end
    end

    %% Serial MPC
    a0=a01;
    x0=x01;
    x_eq=x_eq1;
    tmeasure=t_0;
    xmeasure=[0;0;0;0];

    %Define subsystem dimensions
    n=2;
    m=1;
    % stage costs - different for each subsystem
    Q1= 0.1*eye(2);
    Q2= 0.01*eye(2);
    %initial and terminal state equality constraints
    Aeq=[eye(n),zeros(n,N*(n+m));zeros(n,N*n) eye(n) zeros(n,N*m)]; %

    %Now enter feedback loop
    for ii = 1:mpciterations % maximal number of iterations

        % Set initial guess and initial constraint
        beq1=[[xmeasure(1);xmeasure(3)];x_eq];
        beq2=[[xmeasure(2);xmeasure(4)];x_eq];

        y_init=[x0;a0];

        t_Start = tic;

        % Solve optimization problem: SUBSYSTEM 1
        %structure: y_OL=[x_OL,u_OL];
        [y_OL, V, exitflag, output]=fmincon(@(y) quadcost( N, y, x_eq, Q1,n,delta),...
            y_init,A,b,Aeq,beq1,lb,ub,...
            @(y) nonlinearconstraints_distrib(N, delta,y,yneigh,n,m), options);

        t_Elapsed1 = toc( t_Start ); 

        x_OL1 = y_OL(1:n*(N+1));
        a_OL1 = y_OL(n*(N+1)+1:end);


        %Perform communication
        yneigh=y_OL; %contains the optimally predicted states and controls for the subsystem that was just optimized

        %Initial guess for subsystem 2 optimization: 
        a0 = [a_OL2(m+1:end); zeros(m,1)];
        x0 = [x_OL2(n+1:end); Decoupleddynamic_distrib(delta, x_OL2(end-n+1:end), a_OL2(end-m+1:end))];

        y_init=[x0;a0];

         t_Start = tic;

        % Solve optimization problem: SUBSYSTEM 2
        %structure: y_OL=[x_OL,u_OL];


        [y_OL, V, exitflag, output]=fmincon(@(y) quadcost( N, y, x_eq, Q2,n,delta),...
            y_init,A,b,Aeq,beq2,lb,ub,...
            @(y) nonlinearconstraints_distrib(N, delta, y,yneigh,n,m), options);

        t_Elapsed2 = toc( t_Start ); 

        x_OL2=y_OL(1:n*(N+1));
        a_OL2=y_OL(n*(N+1)+1:end);

        yneigh=y_OL;

        %Initial guess for subsystem 1 optimization: 
        a0 = [a_OL1(m+1:end); zeros(m,1)];
        x0 = [x_OL1(n+1:end); Decoupleddynamic_distrib(delta, x_OL1(end-n+1:end), a_OL1(end-m+1:end))];

        %Calculate actual torque inputs
        taumeasure=computeTau(xmeasure,[a_OL1(1);a_OL2(1)]); 

        % Store closed loop data
        t = [ t, tmeasure ];
        x = [ x, xmeasure ];
        a = [ a, [a_OL1(1);a_OL2(1)] ];
        tau=[tau,taumeasure];

        %Update centralized system with the 2 controls
        xmeasure = Decoupleddynamic(delta,xmeasure,[a_OL1(1);a_OL2(1)]);
        tmeasure = tmeasure + delta;   

        %
        % Print numbers
        fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f  %+6.3f %+6.3f\n', ii, tau(1,end),...
                tau(2,end),x(1,end),x(2,end),t_Elapsed1,t_Elapsed2);

        %plot closed-loop state and control trajetories    
        %figure(1);
        subplot(2,2,1)
        plot(tmeasure,xmeasure(1),'b*')
        subplot(2,2,2)
        plot(tmeasure,taumeasure(1),'b*');hold on
        subplot(2,2,3)
        %figure(2);
        plot(tmeasure,xmeasure(2),'r*')
        subplot(2,2,4)
        plot(tmeasure,taumeasure(2),'r*');hold on
        drawnow

    end
    subplot(2,2,1)
    title('Serial DMPC: 2-link robot joints','Fontsize',14)
    ylabel('q_1','Fontsize',14)
    legend('Centralized','Serial DMPC')
    subplot(2,2,3)
    ylabel('q_2','Fontsize',14)
    xlabel('Time,sec','Fontsize',14)
    subplot(2,2,2)
    ylabel('\tau_1')
    title('Serial DMPC: 2-link robot torques','Fontsize',14)
    subplot(2,2,4)
    ylabel('\tau_2','Fontsize',14)
    xlabel('Time,sec','Fontsize',14)
end  
  

function [yp,yneigh]=splitdata(y,p,N)
    %This function extracts neighbor information from a centralized feasible
    %solution
    %Works only for 2 subsystems each with 2 states and 1 control! 
    yp=zeros(2*(N+1)+N,1);
    yneigh=zeros(2*(N+1)+N,1);
    if p==1,
        %extract odd entries for yp and even for yneigh 
        for i=0:2*N+1+N,
            yp(i+1)=y(2*i+1);
            yneigh(i+1)=y(2*(i+1));
        end
    elseif p==2
         %extract odd entries for yneigh and even for yp
        for i=0:2N+1+N,
            yneigh(i+1)=y(2*i+1);
            yp(i+1)=y(2*(i+1));
        end
    else
        error('Bad subsystem index')
    end
end
          

function [x] = dynamic(delta, x0, u)
    %discretize the system using a Runge-Kutta scheme (ode45)
    
    % Set options
    atol_ode  = 1e-4;
    rtol_ode  = 1e-4;
    options = odeset('AbsTol', atol_ode, 'RelTol', rtol_ode);
    
    % Evaluate the system dynamics (integration)
    [t_intermediate, x_intermediate] = ode23(@system, [0,delta], x0, options, u);
        
    x = x_intermediate(end,:)';

end

function [x] = Decoupleddynamic(delta, x0, a)
    %discretize the system using a Runge-Kutta scheme (ode45)
    
    % Set options
    atol_ode  = 1e-4;
    rtol_ode  = 1e-4;
    options = odeset('AbsTol', atol_ode, 'RelTol', rtol_ode);
    
    % Evaluate the system dynamics (integration)
    [t_intermediate, x_intermediate] = ode23(@Decoupledsystem, [0,delta], x0, options, a);
        
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


function u = computeTau(x, a)
    % System dynamics: two-link planar manipulator in state-space

    q1=x(1);q1dot=x(3);q2=x(2);q2dot=x(4); %parse
    u = zeros(2,1);
    
    %Physical parameters
    
    I1=2;
    I2=1;
    lc1=0.5;
    l1=1;
    lc2=0.5;
    m1=1;
    m2=0.5;
    g=9.8;
    %System matrices
    
    M=[ I1 + I2 + lc1^2*m1 + m2*(l1^2 + lc2^2) + 2*l1*lc2*m2*cos(q2), m2*lc2^2 + l1*m2*cos(q2)*lc2 + I2;
                           m2*lc2^2 + l1*m2*cos(q2)*lc2 + I2,                     m2*lc2^2 + I2];    
    C=[ -l1*lc2*m2*q2dot*sin(q2), - l1*lc2*m2*q1dot*sin(q2) - l1*lc2*m2*q2dot*sin(q2);
    l1*lc2*m2*q1dot*sin(q2),                                                   0];
    gv= [g*cos(q1)*(l1*m2 + lc1*m1) + g*lc2*m2*cos(q1 + q2);
                              g*lc2*m2*cos(q1 + q2)];
   
    u=M*a+C*[q1dot;q2dot]+gv;
end

function xdot = Decoupledsystem(~, x, a)
    % System dynamics for decoupled double integrators
    % Used in optimization
    q1=x(1);q1dot=x(3);q2=x(2);q2dot=x(4); %parse
    xdot = zeros(4,1); %pre-allocate
    xdot=[x(3:4);a];    
end

function xdot = Decoupledsystem_distrib(~, x, a)
    % System dynamics for a double integrator
    % Used in optimization
    q1=x(1);q1dot=x(2);
    xdot=[x(2);a];    
end



function cost = quadcost(N, y, x_eq, Q, n,delta)
%NOTE: this same function can be used in the distributed case if the
%passed arguments are suitable adjusted (for instance, y contains only the
%states and control of the given subsystem, Q is a block diagonal component
%of the "total" Q, etc.

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
   a=y(n*(N+1)+1:end);
   c = [];
   ceq = [];
   % constraints along prediction horizon
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);        
        a_k=a((k-1)*m+1:k*m);
        %dynamic constraint
        ceqnew=x_new - Decoupleddynamic(delta, x_k, a_k);
        ceq = [ceq ceqnew];
        %nonlinear constraints on input : |u1,2|<=UMAX
        %Calculate actual joint torques
        u_k=computeTau(x_k,a_k);    
        %Impose inequality constraints
        UMAX1=80;UMAX2=35; %lower than 35 for UMAX2 may result in infeasibility
        cnew=[u_k;-u_k]-[UMAX1;UMAX2;UMAX1;UMAX2];
        c=[c cnew];
    end
end

function [c, ceq] = nonlinearconstraints_distrib(N, delta, yp, yneigh, n,m) 
   %Evaluates the nonlinear constraints for a subsystem using neighbor
   %information
   %n is the subsystem dimension (assumed constant across subsystems)
   %m is the number of controls of each subsystem (assumed constant across
   %subsystems)
   % Introduce the nonlinear constraints
   x=yp(1:n*(N+1));
   a=yp(n*(N+1)+1:end);
   %*****In this version, there's only 1 neighbor which has the same n and
   %m as its neighbor.
   xneigh=yneigh(1:n*(N+1));
   aneigh=yneigh(n*(N+1)+1:end);
   %*****
   c = [];
   ceq = zeros(2,N);
   % constraints along prediction horizon
    for k=1:N
        x_k=x((k-1)*n+1:k*n);
        x_new=x(k*n+1:(k+1)*n);        
        a_k=a((k-1)*m+1:k*m);
        %dynamic constraint
        ceqnew=x_new - Decoupleddynamic_distrib(delta, x_k, a_k);
        ceq = [ceq ceqnew];
        %nonlinear constraints on input : |u1,2|<=UMAX
        %Calculate actual joint torques -- THIS REQUIRES yneigh
        xneigh_k=xneigh((k-1)*n+1:k*n);
        aneigh_k=aneigh((k-1)*m+1:k*m);
        %Assemble
        xall=[x_k;xneigh_k];
        aall=[a_k;aneigh_k];
        u_k=computeTau(xall,aall);    
        %Impose inequality constraints
        UMAX1=80;UMAX2=35; %lower than 35 for UMAX2 may result in infeasibility
        cnew=[u_k;-u_k]-[UMAX1;UMAX2;UMAX1;UMAX2];
        c=[c cnew];
    end
end
