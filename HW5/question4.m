
% Continuos Plant 
Ac = [0 1; 0 0];
Bc = [0; 1];
Cc = eye(2);
Dc = 0;

% Parameters
Q = eye(2);
R = 1;

delta = 0.1;
T = 1.5;

K0 = place(Ac, Bc, [-1 -2])
x0 = K0;

%%
nonlcon = @(X) nonlconstraints (X, Ac, Bc, Q, R);

options = struct('MaxFunctionEvaluations', 1000, 'MaxIterations', 1000);

fun = @(X) objFunction (X, Ac, Bc, Q, R);
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
[X,FVAL,EXITFLAG,OUTPUT] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

K = X;
Qst = Q + K'*R*K;
Acl = Ac-Bc*K;
P = lyap(Acl', Qst);
Pplot = inv(P);

figure; hold on
E = ellipsoid([0; 0], Pplot);  plot(E)

[X,FVAL,EXITFLAG,OUTPUT] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

hw5_CTdoubleIntMPCregion(P)

function out = objFunction (X, A, B, Q, R)
    K = X;
    Acl = A-B*K;
    Qst = Q + K'*R*K;
    P = lyap(Acl', Qst);
    out = trace(P);
end

function [c, ceq] = nonlconstraints (X, A, B, Q, R)
    c = [];
    K = X;
    Acl = A-B*K;
    Qst = Q + K'*R*K;
    P = lyap(Acl', Qst);
        
    ceq = P(1,2)-P(2,1);
    c = [c; -P(1,1); -P(1,1)*P(2,2)+P(2,1)*P(1,2)];   
end

