
clear all
close all

% getStates();

% function getStates()

    % Load Symbolic Matlab Code
    DE_Generator
    syms u1 u2 u3 u4

    %% Governing Equations
    z1 = [q1; q2; q3; q4];
    z2 = [q1dot; q2dot; q3dot; q4dot];
    z = [z1; z2];
    
    u = [u1; u2; u3; u4];
    
    z1d = z2;
    z2d = inv(M)*(-C*z2 + u);
    zDot = [z1d; z2d];
    
    for i=1:length(z)
        for j=1:length(z)
            DziDzj(i,j) = diff(zDot(i), z(j));
        end
    end
    
    for i=1:length(z)
        for j=1:length(u)
            DziDuj(i,j) = diff(zDot(i), u(j));
        end
    end
    
    %%
    
    z0 = [0; pi/4; 0; pi/4; 0; 0; 0; 0];    zb = z0;
    u0 = [0; 0; 0; 0];                      ub = u0;
    u1 = u0(1);
    u2 = u0(2);
    u3 = u0(3);
    u4 = u0(4);
    
    q1 = z0(1); 
    q2 = z0(2);
    q3 = z0(3);
    q4 = z0(4);
    q1dot = z0(5);
    q2dot = z0(6);
    q3dot = z0(7);
    q4dot = z0(8);
    
    A = double(subs(DziDzj));
    B = double(subs(DziDuj));
    
    Ts = 0.01;
    [Ak, Bk] = c2d(A,B,Ts);
    
     %% Plots
     
     tf = 1;
     sim('openloop.slx', tf);
     
     %%
     figure('Name','Accuracy of the Linearization and the Discretization','units','normalized','outerposition',[0 0 1 1]); 
     subplot(421); hold on; plot(t, q1, 'LineWidth',2);     plot(t, q1lin, 'LineWidth',2);    stairs(dt, q1d, '-k');   xlim([0 tf]); ylabel('q1'); title('Angle Joint 1'); legend('Nonlinear System', 'Linearized System', 'Discrete System');
     subplot(422); hold on; plot(t, q1dot, 'LineWidth',2);  plot(t, q1dotlin, 'LineWidth',2); stairs(dt, q1dotd, '-k');xlim([0 tf]); ylabel('q1dot'); title('Angular Velocity Joint 1'); 
     subplot(423); hold on; plot(t, q2, 'LineWidth',2);     plot(t, q2lin, 'LineWidth',2);    stairs(dt, q2d, '-k');   xlim([0 tf]); ylabel('q2'); title('Angle Joint 2');
     subplot(424); hold on; plot(t, q2dot, 'LineWidth',2);  plot(t, q2dotlin, 'LineWidth',2); stairs(dt, q2dotd, '-k');xlim([0 tf]); ylabel('q2dot'); title('Angular Velocity Joint 2');
     subplot(425); hold on; plot(t, q3, 'LineWidth',2);     plot(t, q3lin, 'LineWidth',2);    stairs(dt, q3d, '-k');   xlim([0 tf]); ylabel('q3'); title('Angle Joint 3');
     subplot(426); hold on; plot(t, q3dot, 'LineWidth',2);  plot(t, q3dotlin, 'LineWidth',2); stairs(dt, q3dotd, '-k');xlim([0 tf]); ylabel('q3dot'); title('Angular Velocity Joint 3');
     subplot(427); hold on; plot(t, q4, 'LineWidth',2);     plot(t, q4lin, 'LineWidth',2);    stairs(dt, q4d, '-k');   xlim([0 tf]); ylabel('q4');     xlabel('Times [s]'); title('Angle Joint 4');
     subplot(428); hold on; plot(t, q4dot, 'LineWidth',2);  plot(t, q4dotlin, 'LineWidth',2); stairs(dt, q4dotd, '-k'); xlim([0 tf]); ylabel('q4dot');  xlabel('Times [s]'); title('Angular Velocity Joint 4');
     
     %% Controller
     
     z0 = zb*1.2;
     
     Q = eye(length(A));
     R = 0.1*eye(length(u));
     K = lqr(A, B, Q, R);
     Kk = dlqr(Ak, Bk, Q, R);
     Kdb = place(Ak,Bk,rand(8,1)*1e-5);
                   
     tf = 10;
     sim('closedloop.slx', tf);
     
     figure('Name','LQR Controller','units','normalized','outerposition',[0 0 1 1]); 
     subplot(421); hold on; plot(t, q1, 'LineWidth',2);     plot(t, q1lin, 'LineWidth',2);    stairs(dt, q1d, '-k');   xlim([0 tf]); ylabel('q1'); title('Angle Joint 1'); legend('Nonlinear System', 'Linearized System', 'Discrete System');
     subplot(422); hold on; plot(t, q1dot, 'LineWidth',2);  plot(t, q1dotlin, 'LineWidth',2); stairs(dt, q1dotd, '-k');xlim([0 tf]); ylabel('q1dot'); title('Angular Velocity Joint 1'); 
     subplot(423); hold on; plot(t, q2, 'LineWidth',2);     plot(t, q2lin, 'LineWidth',2);    stairs(dt, q2d, '-k');   xlim([0 tf]); ylabel('q2'); title('Angle Joint 2');
     subplot(424); hold on; plot(t, q2dot, 'LineWidth',2);  plot(t, q2dotlin, 'LineWidth',2); stairs(dt, q2dotd, '-k');xlim([0 tf]); ylabel('q2dot'); title('Angular Velocity Joint 2');
     subplot(425); hold on; plot(t, q3, 'LineWidth',2);     plot(t, q3lin, 'LineWidth',2);    stairs(dt, q3d, '-k');   xlim([0 tf]); ylabel('q3'); title('Angle Joint 3');
     subplot(426); hold on; plot(t, q3dot, 'LineWidth',2);  plot(t, q3dotlin, 'LineWidth',2); stairs(dt, q3dotd, '-k');xlim([0 tf]); ylabel('q3dot'); title('Angular Velocity Joint 3');
     subplot(427); hold on; plot(t, q4, 'LineWidth',2);     plot(t, q4lin, 'LineWidth',2);    stairs(dt, q4d, '-k');   xlim([0 tf]); ylabel('q4');     xlabel('Times [s]'); title('Angle Joint 4');
     subplot(428); hold on; plot(t, q4dot, 'LineWidth',2);  plot(t, q4dotlin, 'LineWidth',2); stairs(dt, q4dotd, '-k'); xlim([0 tf]); ylabel('q4dot');  xlabel('Times [s]'); title('Angular Velocity Joint 4');
     
     tf = 0.5;
     sim('closedloop.slx', tf);
     
     figure('Name','Deadbeat Controller','units','normalized','outerposition',[0 0 1 1]); 
     subplot(421); hold on; stairs(t, q1d1, 'LineWidth',2);     xlim([0 tf]); ylabel('q1');     title('Angle Joint 1'); 
     subplot(422); hold on; stairs(t, q1dotd1, 'LineWidth',2);	xlim([0 tf]); ylabel('q1dot');  title('Angular Velocity Joint 1'); 
     subplot(423); hold on; stairs(t, q2d1, 'LineWidth',2);     xlim([0 tf]); ylabel('q2');     title('Angle Joint 2');
     subplot(424); hold on; stairs(t, q2dotd1, 'LineWidth',2);	xlim([0 tf]); ylabel('q2dot');  title('Angular Velocity Joint 2');
     subplot(425); hold on; stairs(t, q3d1, 'LineWidth',2);     xlim([0 tf]); ylabel('q3');     title('Angle Joint 3');
     subplot(426); hold on; stairs(t, q3dotd1, 'LineWidth',2);  xlim([0 tf]); ylabel('q3dot');  title('Angular Velocity Joint 3');
     subplot(427); hold on; stairs(t, q4d1, 'LineWidth',2);     xlim([0 tf]); ylabel('q4');     xlabel('Times [s]'); title('Angle Joint 4');
     subplot(428); hold on; stairs(t, q4dotd1, 'LineWidth',2);  xlim([0 tf]); ylabel('q4dot');  xlabel('Times [s]'); title('Angular Velocity Joint 4');
     
% end