
function FiniteHorizon_LQR(A, B, X0, Q, R, Qf, dt, tf)

% Simulation Parametes
N = tf/dt+2;
t = 0:dt:tf;

[Ad, Bd] = c2d(A,B, dt);

% Initialized parameters
P(:,:,1)= Qf;

N = tf/dt+2;
for k=2:N-1
    S = R + Bd'*P(:,:,k-1)*Bd;
    F(:,:,N-k) =  -( inv(S) * Bd' * P(:,:,k-1) * Ad );
    P(:,:,k) = (Ad + Bd*F(:,:,N-k))'*P(:,:,k-1)*(Ad + Bd*F(:,:,N-k)) + F(:,:,N-k)'*R*F(:,:,N-k) + Q;
end

x(1,:) = X0;
for k=1:N-2
    u(k,:) = F(:,:,k)*x(k,:).';    
    XD = A*x(k,:).' + B*u(k,:)';
    
    x(k+1,:) = x(k,:).' + XD*dt;
end

figure
ax1 = subplot(321); plot(ax1, t, x(:,1),'LineWidth',2); 
ax2 = subplot(323); plot(ax2, t, x(:,2),'LineWidth',2);
ax3 = subplot(325); plot(ax3, t, x(:,3),'LineWidth',2); xlabel('Time [s]', 'Interpreter','Latex', 'FontSize',14);
ax4 = subplot(3,2,[2 4 6]); plot3(x(:,1), x(:,3), x(:,3), 'LineWidth',2); axis equal
title(ax1, 'State $x_1$','Interpreter','latex','FontSize',14); 
title(ax2, 'State $x_2$','Interpreter','latex','FontSize',14);   
title(ax3, 'State $x_3$','Interpreter','latex','FontSize',14);   

end
