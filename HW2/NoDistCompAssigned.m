% Erivelton Gualter dos Santos
% 09/22/2018
%
% Homework 2
% Based on the code provided by Hanz Richter at ESC 794 class.

% Disturbance compensation
% clear;
% clc;
% close all

Bd = [1;0];
Cd = eye(2); % For constraints
Cdy = [1 0]; % Output matrix for cost function

% Initial condition
xnow = [1;1];
xnext = xnow;
ynext = Cdy*xnow;
ypred = ynext;

% Horizons 
if ~ exist('guihw')
    ny = 5;
    nu = 4;
    simhor = 8;
    lambda = 2; % Weight   
    UMAX = 0.5; 
end

u = [];
for k=2:simhor+1
    % Update linearization matrices
    x1 = xnow(1); x2 = xnow(2);
    Ad = [cos(x1)*x2 sin(x1); -sin(x1)*x1+cos(x1) 0];% Jacobian at x1, x2
       
    % Build Prediction Matrices for output 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = Cdy; D = 0;
    n = size(Ad,1); % state dimension
    m = size(Bd,2); % input dimension
    p = size(C,1);  % output dimension

    % Compute Py
    P = [];
    P = C*Ad; 
    for i=1:ny-1
      P = [P; C*Ad^(i+1)];
    end
    Py = P;
    
    % Compute Hy
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
    Hy = H;

    % Matrices for unit box constraints on states
    C = Cd; D = [0;0];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pc = size(C,1); % output dimension
    p = pc;
    %Compute Pc
    P = [];
    P = C*Ad; 
    for i=1:ny-1
      P = [P; C*Ad^(i+1)];
    end
    Pc=P;
    % Compute Hc
    H=[];
    H=zeros(p*ny,m*nu);
    for i=1:ny,
        for j=1:i,
            H(1+(i-1)*p:i*p,1+(j-1)*m:j*m) = C*Ad^(i-j)*Bd;
        end
        H(1+(i-1)*p:i*p,1+j*m:(j+1)*m) = D;
    end
    % Retain only the first nu blocks
    H = H(:,1:nu*m);
    Hc = H;

    % Output reference
    r = zeros(ny,1);
    
    %Form the L vectors for output constraints
    Ly=repmat(eye(pc),[ny 1]);

    Ybar=[1;1];
    y_max=Ly*Ybar;
    y_min=-y_max;

    %Control input limits
    UB=UMAX*ones(nu,1);
    LB=-UB;

    S = Hy'*Hy+lambda*eye(m*nu); %Hessian portion

    x=xnow;

    f=2*Hy'*(Py*x-r);

    dyp=y_max-Pc*x;
    dym=-y_min+Pc*x;
    
    M=[Hc;-Hc];
    d=[dyp;dym];

    % Solve quadratic program
    ctrlvec = quadprog(2*S,f,M,d,[],[],LB,UB);
    sol = 1;
    if isempty(ctrlvec)
        sol = 0;
        break
    end
    u_apply = ctrlvec(1:m); % Extract first term of optimal sequence
    ypred_end = Py*xnow+Hy*ctrlvec;
    ypred_end = ypred_end(end);
    ypred = [ypred ypred_end];
    u = [u u_apply];                    % Store control history
    xnext(:,k) = Ad*xnow+Bd*u(:,end);   % Update plant
    ynext = [ynext Cdy*xnext(:,k)];     % Store output history
    xnow = xnext(:,k); 
end

if sol
    ax1 = subplot(221); hold on; plot(xnext(1,:),xnext(2,:)); 
        xlabel('$X_1$', 'interpreter', 'latex', 'FontSize',14);
        ylabel('$X_2$', 'interpreter', 'latex', 'FontSize', 14);

    % figure(2); 
    ax2 = subplot(222); hold on;
    plot(UMAX*ones(size(u)), '-k','LineWidth', 2);  % Constraints
    plot(-UMAX*ones(size(u)), '-k','LineWidth', 2); % Constraints
    plot(u); 
        xlabel('Simulation Horizon', 'interpreter', 'latex', 'FontSize',14);
        ylabel('Control Input', 'interpreter', 'latex', 'FontSize', 14);
        legend('Constraints');

    % figure(3); 
    ax3 = subplot(223); hold on;
    plot(ynext); plot(ypred,'r')
        xlabel('Simulation Horizon', 'interpreter', 'latex', 'FontSize',14);
        ylabel('Output', 'interpreter', 'latex', 'FontSize', 14);
        legend('Output History','Prediction');
end    