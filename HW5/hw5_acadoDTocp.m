clear all;
clc 

BEGIN_ACADO;                               

    acadoSet('problemname', 'hw5_acadoDTocp'); 
    
    DifferentialState       x1 x2;  % The differential states
    Control                 u;      % The controls
       
    %% Differential Equation
    f = acado.DiscretizedDifferentialEquation(1);
    
    f.add(next(x1) == x1*(1-u));    
    f.add(next(x2) == sqrt(x1*x1+x2*x2)*u);
    
    %% Optimal Control Problem
    ocp = acado.OCP([0:2]);
                                            
    ocp.minimizeLagrangeTerm(x1*x1+x2*x2);   % Minimize a Lagrange term
    
    ocp.subjectTo( f );     
    
    ocp.subjectTo( 'AT_START', x1 == 1.0 );  % Initial condition
    ocp.subjectTo( 'AT_START', x2 == 1.0 ); 
    ocp.subjectTo(  0 <= u <= 1 );
    ocp.subjectTo( 'AT_END', x1 == 0);
    ocp.subjectTo( 'AT_END', x2 == 0);
     
    %% Optimization Algorithm
    algo = acado.OptimizationAlgorithm(ocp); % Set up the optimization algorithm
    
    algo.set('INTEGRATOR_TOLERANCE', 1e-5 ); % Set some parameters for the algorithm
   
END_ACADO;           

out = hw5_acadoDTocp_RUN();         
                                            
%% Plots
figure('Name','Question 5','units','normalized', ...
    'outerposition',[0 0 1 1], 'NumberTitle','off');
subplot(221); stairs(out.STATES(:,1), out.STATES(:,2), 'LineWidth', 2); xlabel('Time'); ylabel('x1'); title('State Response x1');
subplot(222); stairs(out.STATES(:,1), out.STATES(:,3), 'LineWidth', 2); xlabel('Time'); ylabel('x2'); title('State Response x2');
% out.STATES(:,1) : Time Points
% out.STATES(:,2) : State x

subplot(212)
stairs(out.STATES(:,1), out.CONTROLS(:,2), 'r', 'LineWidth', 2); xlabel('Time'); ylabel('u'); title('Control Input');

u1=out.CONTROLS(:,2); %save for comparisons below
