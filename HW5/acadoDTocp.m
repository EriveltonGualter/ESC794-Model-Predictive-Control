%clear;

BEGIN_ACADO;                                % Always start with "BEGIN_ACADO". 
    
    acadoSet('problemname', 'acadoDTocp');     % Set your problemname. If you 
                                                    % skip this, all files will
                                                    % be named "myAcadoProblem"
    
    
    DifferentialState       x1 x2;              % The differential states
    Control                 u;              % The controls
    %Disturbance             w;              % The disturbances
    %Parameter               p q;            % The free parameters
    
    
    
    %% Differential Equation
    f = acado.DiscretizedDifferentialEquation(1);       % Set the differential equation object
    
    f.add(next(x1) == x1+2*x2);    
    f.add(next(x2) == -x2+2*u);
                                            % Write down your ODE. You can add 
                                            % multiple ODE's with f.add() when you 
                                            % have multiple states. To print an ODE 
                                            % to the screen, use
                                            % f.differentialList{1}.toString

    
    %% Optimal Control Problem
    ocp = acado.OCP([0:3]);          % Set up the Optimal Control Problem (OCP)
                                            % Start at 0s, control in 3
                                            % intervals up to 3s
                                            
    ocp.minimizeLagrangeTerm(u^2+(x1-4)^2+x2^2);   % Minimize a Lagrange term
    
    ocp.subjectTo( f );                     % Your OCP is always subject to your 
                                            % differential (or difference) equation
                                            
    ocp.subjectTo( 'AT_START', x1 == 0.0 );  % Initial condition
    ocp.subjectTo( 'AT_START', x2 == 0.0 ); 
    %ocp.subjectTo(  0.1 <= p <= 2.0 );      % Bounds
    %ocp.subjectTo(  -1 <= u <= 1 );
    %ocp.subjectTo( -0.1 <= w <= 2.1 );
    %ocp.subjectTo('AT_END',x1==4);
    %ocp.subjectTo('AT_END',x2==0);
    
    %% Optimization Algorithm
    algo = acado.OptimizationAlgorithm(ocp); % Set up the optimization algorithm

    
    algo.set('INTEGRATOR_TOLERANCE', 1e-5 ); % Set some parameters for the algorithm
    
    
END_ACADO;           % Always end with "END_ACADO".
                     % This will generate a file problemname_ACADO.m. 
                     % Run this file to get your results. You can
                     % run the file problemname_ACADO.m as many
                     % times as you want without having to compile again.



out = acadoDTocp_RUN();                % Run the test. The name of the RUN file
                                            % is problemname_RUN, so in
                                            % this case getting_started_RUN
                                            
%draw
subplot(2,1,1)
stairs(out.STATES(:,1), out.STATES(:,2), 'r')
% The first column contains the time points
% The second column contains the state 'x'
title('x');

subplot(2,1,2)
stairs(out.STATES(:,1), out.CONTROLS(:,2), 'r')
title('u');

u1=out.CONTROLS(:,2); %save for comparisons below

%Now illustrate the Dynamic Programming Principle:
%Solve again with the state at k=1 as the initial state, and N=2

%The previous optimal cost and stage 1 cost were 1.5 and 0.25

BEGIN_ACADO;                                % Always start with "BEGIN_ACADO". 
    
    acadoSet('problemname', 'acadoDTocp');     % Set your problemname. If you 
                                                    % skip this, all files will
                                                    % be named "myAcadoProblem"
    
    
    DifferentialState       x1 x2;              % The differential states
    Control                 u;              % The controls
    %Disturbance             w;              % The disturbances
    %Parameter               p q;            % The free parameters
    
    
    
    %% Differential Equation
    f = acado.DiscretizedDifferentialEquation(1);       % Set the differential equation object
    
    f.add(next(x1) == x1+2*x2);    
    f.add(next(x2) == -x2+2*u);
                                            % Write down your ODE. You can add 
                                            % multiple ODE's with f.add() when you 
                                            % have multiple states. To print an ODE 
                                            % to the screen, use
                                            % f.differentialList{1}.toString

    
    %% Optimal Control Problem
    ocp = acado.OCP([0:2]);          % Set up the Optimal Control Problem (OCP)
                                            % Start at 0s, control in 21
                                            % intervals up to 2s
                                            
    ocp.minimizeLagrangeTerm(u^2);   % Minimize a Lagrange term
    
    ocp.subjectTo( f );                     % Your OCP is always subject to your 
                                            % differential (or difference) equation
                                            
    ocp.subjectTo( 'AT_START', x1 == 0.0 );  % Initial condition
    ocp.subjectTo( 'AT_START', x2 == 1.0 ); 
    %ocp.subjectTo(  0.1 <= p <= 2.0 );      % Bounds
    %ocp.subjectTo(  -1 <= u <= 1 );
    %ocp.subjectTo( -0.1 <= w <= 2.1 );
    %ocp.subjectTo('AT_END',x1==4);
    %ocp.subjectTo('AT_END',x2==0);
    
    %% Optimization Algorithm
    algo = acado.OptimizationAlgorithm(ocp); % Set up the optimization algorithm

    
    algo.set('INTEGRATOR_TOLERANCE', 1e-5 ); % Set some parameters for the algorithm
    
    
END_ACADO;

tic
out = acadoDTocp_RUN(); 
toc
%The optimal cost-to-go is 1.25
%Check 1.5=0.25+1.25

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also check that this solution is just the tail of the previous one
u_tail=out.CONTROLS(:,2)
u1

