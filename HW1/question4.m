
clear all
close all

% Continuos Plant 
A = [0 1; 0 0];
B = [0; 1];
C = eye(2);
D = 0;

% Discrete Plant
Ts = 0.01;
sys = ss(A,B,C,D);
sysd = c2d(sys,Ts,'zoh');
Ad = sysd.a;
Bd = sysd.b;
Cd = sysd.c;
Dd = sysd.d;

[Ad, Bd, Cd, Dd] = ssdata(sysd);

p = [0 0];
K = dlqr(Ad, Bd, diag([100 1]),0.1);
Acl = Ad-Bd*K;

Q = eye(2);

P = dlyap(Acl', Q');
%%
nonlcon = @(X) nonlconstraints (X, Acl);
x0 = reshape(Q,4,1);

options = struct('MaxFunctionEvaluations', 1000, 'MaxIterations', 1000);

fun = @(X) objFunction (X, Acl);
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

Qnew = reshape(x,2,2);
Pnew = dlyap(Acl', Qnew');
figure; hold on
rectangle('Position', [-1 -1 2 2]); 
Pnew = inv(Pnew);
E = ellipsoid([0; 0], Pnew); 
plot(E)
axis equal

[x10,x20] = ginput(1);
X = [x10; x20];

t = 0:Ts:10;
sys_cl = ss(Ad-Bd*K,B,C,D,Ts);
[X, t] = lsim(sys_cl,zeros(size(t)), t, X);

plot(X(:,1), X(:,2))

function out = objFunction (X, Acl)
    Q = reshape(X, 2,2);
    P = dlyap(Acl', Q');
    out = trace(P);
end

function [c, ceq] = nonlconstraints (X, Acl)
    c = [];
    Q = reshape(X, 2,2);

    P = dlyap(Acl', Q');
    gama1 = sqrt([1 0]*inv(P)*[1 0]');
    gama3 = sqrt([0 1]*inv(P)*[0 1]');
        
    ceq = X(2)-X(3);
    c = [c; -Q(1,1); -Q(1,1)*Q(2,2)+Q(2,1)*Q(1,2); gama1-1; gama3-1];   
end

