clear all
close all

% System parameters
Jt  = 20;   % total inertia reflected to linear coordinate, kg-m^2
a   = 12;   % gear ratio times torque constant, N-m/A
Rm  = 0.1;  % motor resistance, Ohm
R   = 1;    % drum radius, m
mc  = 1000; % nominal elevator car mass, with passengers, kg
mw  = 750;  % counterweight mass, kg
g   = 9.81; % gravity, m/s^2

mcw  = (mc - mw);

% Continuous-time State-Space matrices with pos - vel states
A = [0 1; 0 -a^2/(Rm*Jt)];  B = [1; 1/Jt];
C = eye(2); D = 0;

tspan = [0 2];
u = 0;
X0 = [1 0];
Ts = 0.05;

sys = ss(A, B, C, D);
tl = tspan(1):Ts:tspan(2);
f = 1;
[Xl, tl] = lsim(sys, sin(2*pi*tl*f), tl);

sysd = c2d(sys, Ts, 'zoh');
[Ad,Bd,Cd,Dd] = ssdata(sysd);
[Xnl, tn] = lsim(sysd, sin(2*pi*tl*f), tl);

% Plots
subplot(211); hold on; stairs(tn, Xnl(:,1), 'LineWidth', 1); plot(tl, Xl(:,1), 'LineWidth', 2); 
    ylabel('$x_1$', 'interpreter', 'latex', 'FontSize', 18);
subplot(212); hold on; stairs(tn, Xnl(:,2), 'LineWidth', 1); plot(tl, Xl(:,2), 'LineWidth', 2); 
    ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 18); 
    xlabel('Time [s]', 'interpreter', 'latex', 'FontSize', 16);