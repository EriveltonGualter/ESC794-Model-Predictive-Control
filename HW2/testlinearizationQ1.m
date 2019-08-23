% First of all lets check linearity 

clear all
close all


tspan = [0 10];
u = 0;
X0 = [1 0];

% Find equilibrium points
ub = 0;
x1b = acos(0);
x2b = 0;
X0 = [x1b; x2b];

% Nonlinear Plant Simulation
f = 1; % Hz
fun = @(t, X) odefun_nl(t, X, sin(2*pi*t*f));
[tn, Xnl] = ode45(fun, tspan, X0);

% Linear Plant Simualation
x1 = x1b; x2 = x2b;
A = [cos(x1)*x2 sin(x1); -sin(x1)*x1+cos(x1) 0];
B = [1; 0];
C = eye(2);
D = 0;
    
sys = ss(A, B, C, D);
tl = linspace(tspan(1),tspan(2), 100);
[Xl, tl] = lsim(sys, sin(2*pi*tl*f), tl);

% Plots
subplot(211); hold on; plot(tn, Xnl(:,1), 'LineWidth', 2); plot(tl, Xl(:,1)+x1b,'--', 'LineWidth', 2); 
    ylabel('$x_1$', 'interpreter', 'latex', 'FontSize', 18);
    title('Test Linearization at Equilibrium point', 'Interpreter', 'latex','FontSize', 16);
subplot(212); hold on; plot(tn, Xnl(:,2), 'LineWidth', 2); plot(tl, Xl(:,2)+x2b,'--', 'LineWidth', 2); 
    ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 18); 
    xlabel('Time [s]', 'interpreter', 'latex', 'FontSize', 16);

function dxdt = odefun_nl(t, X, u)
    x1 = X(1);
    x2 = X(2);
    
    xd1 = sin(x1)*x2 + u;
    xd2 = cos(x1)*x1;
    
    dxdt = [xd1; xd2];
end