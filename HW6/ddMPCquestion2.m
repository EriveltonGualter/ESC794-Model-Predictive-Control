% Parallel distributed MPC

% Erivelton Gualter 
%
% Based on ddMPCexampleUconst.m by Dr. Hanz Richter

function ddMPCquestion2

    clear; clc;
    
    %% Centralized MPC
    % Decoupled plant
    Ad = [1 2 0;0 -1 0;0 0 -2];
    Bd = [1 0;-1 0;0 1];

    %Objective: regulation to an equilibrium state
    %Find some eq state for a specific equilibrium control
    ueq = [0;0];
    xeq = [0; 0; 0];
    
    % Inicial Condition
    X0 = [-.3; .3; .3];

    n = size(Ad,1); %state dimension
    p = n;
    m = size(Bd,2); %input dimension

    % Horizon
    ny = 5;
    nu = ny;

    %Solve and simulate centralized problem
    Q = eye(3);
    R = eye(2);

    [Px, Hx] = buildPH(Ad,Bd,nu,ny,n);

    %Set up matrices for constraints on control

    Cconstr = eye(3); %dummy
    pc = size(Cconstr,1); %dummy
    [Pc, Hc, Cc] = buildPHC(Cconstr,Ad,Bd,ny,nu,pc);

    xeq_hat = xeq;
    ueq_hat = ueq;
    Qbar = Q;
    Rbar = R;

    u_max = [1; 1];
    u_min = [-1; -1];

    [M,d,Qbar,Rbar,xeq_hat,ueq_hat]=buildMdEq(Cc,xeq,ueq,Q,R,u_max,u_min,ny,nu,pc);

    %Form cost function - constant portion
    S = 2*(Hx'*Qbar*Hx+Rbar);
    simhor = 10; 

    %Simulation (unit step reference)
    %Initial conditions
    xa = X0;

    %initialize variable part of cost
    f = 2*((xa'*Px'-xeq_hat')*Qbar*Hx-ueq_hat'*Rbar);

    uprev = zeros(m,1);
    xnow = zeros(3,1);
    u = uprev;

    % Set up equality constr matrices
    Qc  = [zeros(n,n*(ny-1)) eye(n)]*Px;
    Aeq = [zeros(n,n*(ny-1)) eye(n)]*Hx;
    %beq defined inside mpc loop

    T=[];
    for k=1:simhor
       beq = xeq-Qc*xa;
       t = tic;
       ctrlvec = quadprog(S,f,M,d,Aeq,beq);
       T(k) = toc(t);
       u_apply = ctrlvec(1:m); %extract first term of optimal sequence
       u = [u u_apply];  %store control history
       xnext(:,k) = Ad*xnow+Bd*u(:,end); %update plant
       xnow = xnext(:,k); 
       uprev=u(:,end);
       xa = xnow; 
       f = 2*((xa'*Px'-xeq_hat')*Qbar*Hx-ueq_hat'*Rbar);
    end

    figure(4);
    subplot(311); stairs(xnext(1,:),'LineWidth',2); ylabel('z1'); legend('Centralized','Distrib'); title('States: Centralized vs. Parallel MPC','Fontsize',12)
    subplot(312); stairs(xnext(2,:),'LineWidth',2); ylabel('z2'); legend('Centralized','Distrib'); 
    subplot(313); stairs(xnext(3,:),'LineWidth',2); xlabel('Number of Iterations'); ylabel('z3'); legend('Centralized','Distrib'); 

    figure(5);
    subplot(211); stairs(u(1,:),'LineWidth',2); xlabel('Number of Iterations'); ylabel('Control Input, u1'); title('Control Inputs'); xlim([1 simhor]);
    subplot(212); stairs(u(2,:),'LineWidth',2); xlabel('Number of Iterations'); ylabel('Control Input, u2'); xlim([1 simhor]);

    %% Distributed MPC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now solve a dual decomposition MPC with the partition {x1 x2}, {x3}

    % Reset plant state
    xa = X0;

    % Sub-problems
    A11 = Ad(1:2,1:2); 
    A22 = Ad(3,3); 
    A12 = Ad(1:2,3);
    A21 = Ad(3,1:2);  

    B1 = Bd(1:2,1); 
    B2 = Bd(3,2); 
    Q1c = Qc(:,1:2);
    Q2c = Qc(:,3);

    R1c = Aeq(:,[1 3 5 7 9]);  % odd cols of Aeq
    R2c = Aeq(:,[2 4 6 8 10]); % even cols of Aeq

    % Prediction matrices for Subsystem 1: (A11, B1) %%%%%%%%%%%%%%%%%%%%%%
    n = size(A11,1);
    m = size(B1,2);
    [Px1, Hx1] = buildPH(A11,B1,nu,ny,n);

    % Matrices for constraints on control
    % Compute Pc, Hc and Cc for constraints
    Cd_constr1 = eye(2); %dummy
    pc1 = size(Cd_constr1,1); %dummy

    [~, ~, Cc1] = buildPHC(Cd_constr1,A11,B1,ny,nu,pc1);

    % Compute M and d
    xeq1 = xeq(1:2);
    ueq1 = ueq(1);
    umax1 = u_max(1);
    umin1 = u_min(1);

    % LQR gains for Sybsystem 1
    R1 = 1;
    Q1 = eye(2); %for local optimization

    [M1,d1,Qbar1,Rbar1,xeq1_hat,ueq1_hat] = ...
        buildMdEq(Cc1,xeq1,ueq1,Q1,R1,umax1,umin1,ny,nu,pc1);
%     legend('Centralized','Distrib')
%     subplot(2,1,1)
%     ylabel('u_1')
%     legend('Centralized','Distrib')

    % Form cost function - quadratic term 
    S1 = 2*Hx1'*Qbar1*Hx1+Rbar1;

    xa1 = xa(1:2);

    % Prediction matrices for Subsystem 2: (A22, B2) %%%%%%%%%%%%%%%%%%%%%%
    n = size(A22,1);
    [Px2, Hx2] = buildPH(A22,B2,nu,ny,n);

    % Set up matrices for constraints on control
    % Compute Pc, Hc and Cc for constraints
    Cd_constr2 = eye(1);
    pc2 = size(Cd_constr2,1);
    
    [~, ~, Cc2] = buildPHC(Cd_constr2, A22, B2, ny, nu, pc2);

    %Compute M and d
    xeq2 = xeq(2);
    ueq2 = ueq(2);
    umax2 = u_max(2);
    umin2 = u_min(2);

    % LQR gains for Sybsystem 2
    R2 = 1;
    Q2 = 1;

    [M2, d2, Qbar2, Rbar2, xeq2_hat, ueq2_hat] = ...
        buildMdEq(Cc2,xeq2,ueq2,Q2,R2,umax2,umin2,ny,nu,pc2);

    % Form cost function - constant portion
    S2 =2 *Hx2'*Qbar2*Hx2+Rbar2;

    xa2 = xa(3);

    %Prepare for MPC loop

    %Initialize global Lagrange multiplier
    Lambda = [0;0;0];
%     Lambda = rand(3,1);
    cr = 2.5; %convergence rate

    m = size(Bd,2);
    uprev = zeros(m,1);
    xnow = zeros(3,1);
    u = uprev;

    %Main MPC loop
    nEC = []
    for k=1:simhor
       T1 = 0; T2 = 0; 
       converged = 0;
       xa1 = xa(1:2);
       xa2 = xa(3);

       while ~converged 

           f1 = 2*((xa1'*Px1'-xeq1_hat')*Qbar1*Hx1-ueq1_hat'*Rbar1)+Lambda'*R1c;
           f2 = 2*((xa2'*Px2'-xeq2_hat')*Qbar2*Hx2-ueq2_hat'*Rbar2)+Lambda'*R2c;

           t1 = tic;
           Y1 =quadprog(S1,f1,M1,d1);
           dt1 = toc(t1);
           T1 = T1+dt1;

           t2 = tic;
           Y2 = quadprog(S2,f2,M2,d2); 
           dt2 = toc(t2);
           T2 = T2+dt2;
 
           %Find constraint error

           %interleave Y1 and Y2 solutions HARD-CODED FOR GIVEN DIMENSIONS AND
           %HORIZON 

           Y = [Y1(1);Y2(1);Y1(2);Y2(2);Y1(3);Y2(3);Y1(4);Y2(4);Y1(5);Y2(5)];
           ec=Qc*xa+Aeq*Y-xeq;

           % Central Lagrangian update
           Lambda=Lambda+cr*ec; 
           normEC=norm(ec)
           converged = normEC<0.001
           nEC = [nEC normEC];
       end

       Ttotal1(k)=T1; %total time spent by Lagrangian iterations, for each "processor"
       Ttotal2(k)=T2;

       u_apply=[Y1(1);Y2(1)]; %extract first term of optimal sequence
       u=[u u_apply];  %store control history
       xnext(:,k)=Ad*xnow+Bd*u(:,end); %update plant
       xnow=xnext(:,k); 
       uprev=u(:,end);
       xa=xnow;
    end

    figure(4); 
    subplot(3,1,1); hold on; stairs(xnext(1,:),'r--','LineWidth',2);
    subplot(3,1,2); hold on; stairs(xnext(2,:),'r--','LineWidth',2); 
    subplot(3,1,3); hold on; stairs(xnext(3,:),'r--','LineWidth',2)
    
    figure(5);
    title('Control Inputs: Centralized vs. Parallel MPC - Input Constraints')
    subplot(2,1,1); hold on; stairs(u(1,:),'r--'); legend('Centralized','Distrib')
    subplot(2,1,2); hold on; stairs(u(2,:),'r--'); legend('Centralized','Distrib')
    
    
    figure(6)
    bar([Ttotal1; Ttotal2; T]');     
    
    legend('Distrib, processor 1','Distrib, processor 2','Centralized quadprog')
    title('Computation Time: Centralized vs. Parallel MPC - Input Constraints')
    ylabel('Time per quadratic optimization')
    xlabel('MPC iteration')
    legend('Distrib, processor 1','Distrib, processor 2','Centralized quadprog')

%     figure
%     plot(nEC)
    end

    function [Px,Hx]=buildPH(A,B,nu,ny,p)

    n=size(A,1);m=size(B,2);

    %Compute Px
    Px=A; 
    for i=1:ny-1
      Px=[Px;A^(i+1)];
    end

    %Compute Hx
    Hx=zeros(p*ny,m*nu);
    for i=1:ny,
        for j=1:i,
            Hx(1+(i-1)*p:i*p,1+(j-1)*m:j*m)=A^(i-j)*B;
        end
    end

end


function [Pc,Hc,Cc]=buildPHC(Cconstr,A,B,ny,nu,pc)

    m=size(B,2);

    Pc=Cconstr*A;
    for i=1:ny-1,
      Pc=[Pc;Cconstr*A^(i+1)];
    end

    %Compute H for constraints

    Hc=zeros(pc*ny,m*nu);
    for i=1:ny,
        for j=  1:i,
            Hc(1+(i-1)*pc:i*pc,1+(j-1)*m:j*m)=Cconstr*A^(i-j)*B;
        end
    end

    %Compute Cc 
    Cc=zeros(m*nu,m*nu);
    for i=1:nu,
        for j=1:i,
            Cc(1+(i-1)*m:i*m,1+(j-1)*m:j*m)=eye(m);
        end
    end

end

function [M,d,Qbar,Rbar,xeq_hat,ueq_hat]=buildMdEq(Cc,xeq,ueq,Q,R,umax,umin,ny,nu,pc)

    Ly = eye(pc);
    xeq_hat = xeq;
    ueq_hat = ueq;
    Qbar = Q;
    Rbar = R;
    m = size(R,1);

    for i=1:ny-1,
        Ly = [Ly;eye(pc)];
        xeq_hat = [xeq_hat;xeq];
        ueq_hat = [ueq_hat;ueq];
        Qbar = blkdiag(Qbar,Q);
        Rbar = blkdiag(Rbar,R);
    end

    Lu = eye(m);
    for i=1:nu-1,
        Lu = [Lu;eye(m)];
    end

    dup = Lu*umax;
    dum = -Lu*umin;

    M = [Cc;-Cc];
    d = [dup;dum];
end


