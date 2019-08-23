% Erivelton Gualter 

clear all
close all
clc

% Ax>=b
A = [1 1 1; 0 0 -1; 0 -1 0; -1 0 0]; 
B = [4; -2; 0; -4];

P = diag([2 4 8]);
xc = [1 1 1]';
f = -xc'*P;
X = quadprog(P,f,-A,-B)

hold on;
plot3(X(1),X(2),X(3),'s')

% E = ellipsoid(xc, P);
% plot(E,'b');
% 
% Form hyperplane array (box constraints)
HA1 = hyperplane(A(1,:)',B(1));
HA2 = hyperplane(-A(2,:)',-B(2));
HA3 = hyperplane(-A(3,:)',-B(3));
HA4 = hyperplane(-A(4,:)',-B(4));

% plot(HA1,'r')
% plot(HA2,'g')
% plot(HA3,'b')
% plot(HA4,'c')

axis([-4 4 -4 4 -4 4])

x0 = [3.5 -0.5 2]';

plot3(x0(1),x0(2),x0(3),'ok')
plot(HA2,'g')
Xr = x0; 

A*x0-B

% Constraints 1
A2=A(2,:);
Z=null(A2);

g = P*x0 + f';
pz=-(Z'*P*Z)\(Z'*g);
p0=Z*pz

xtest=x0+p0;
A*xtest-B

alph_max3=(B(3)-A(3,:)*x0)/(A(3,:)*p0);
alph_max=alph_max3;
x1=x0+alph_max*p0;

Xr = [Xr x1]; 
plot3(Xr(1,:), Xr(2,:), Xr(3,:),'LineWidth',3)
plot(HA3,'b')
plot3(x1(1),x1(2),x1(3),'*k')

g = P*x1 + f'

% Constraints 2
A3=A(3,:);
Z=null(A3);

g = P*x1 + f';
pz=-(Z'*P*Z)\(Z'*g);
p1=Z*pz

xtest=x1+p1;
A*xtest-B

alph_max1=(B(1)-A(1,:)*x1)/(A(1,:)*p1);
alph_max=alph_max1;
x2=x1+alph_max*p1;

plot(HA1,'r')
plot3(x2(1),x2(2),x2(3),'*k')
Xr = [Xr x2];
plot3(Xr(1,:), Xr(2,:), Xr(3,:),'LineWidth',3)

% Constraints 3
A14=A([1 4],:);
Z=null(A14);

g = P*x2 + f';
pz=-(Z'*P*Z)\(Z'*g);
p2=Z*pz

xtest=x2+p2;
A*xtest-B

alph_max4=(B(4)-A(4,:)*x2)/(A(4,:)*p2);
alph_max1=(B(1)-A(1,:)*x2)/(A(1,:)*p2);
alph_max3=(B(3)-A(3,:)*x2)/(A(3,:)*p2);
alph_max=min(alph_max1,alph_max3);
x3=x2+alph_max*p2;

H = P;
g = f;
lb = [];
ub = [];
lbA = [];
ubA = B;

% [x,fval,exitflag,iter,lambda,auxOutput] = qpOASES( H,g',A,lb,ub,lbA,ubA)

X = quadprog(P,-xc'*P,-A,-B)
[x,fval] = qpOASES( H,g',-A,lb,ub,[],-ubA )

%%
% qpOASES( H,g,A,lb,ub,lbA,ubA{,options{,auxInput}} )
% 
%  min 0.5*x'*H*x + f'*x   subject to:  A*x <= b 
% 
%  min   0.5*x'Hx + x'g
%                  s.t.  lb  <=  x <= ub
%                        lbA <= Ax <= ubA  {optional}
%  
 