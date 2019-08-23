A=[-1 1;-2 -3;-4 3;2 1]; B=[2;-6;-12;-4];
%(In the Ax>=b format)

P=[1 2;2 6];
%Obj fcn is (1/2)x'*P*x

X=quadprog(P,[0 0],-A,-B)

%Draw
E=ellipsoid(inv(P));
% E=ellipsoid(0.7273*inv(P));
plot(E,'b'); hold on
%Form hyperplane array (box constraints)
HA1=hyperplane([-1;1],2);
HA2=hyperplane([2;3],6);
HA3=hyperplane([4;-3],12);
HA4=hyperplane([2;1], -4);

plot(HA1,'r')
plot(HA2,'b')
plot(HA3,'g')
plot(HA4,'k')

axis([-4 4 -4 4])
axis equal

%Attempt a manual active set solution

%Guess an active set: constraint 2
%

A2=A(2,:);
%Find the null space basis
Z=null(A2);
%Select a feasible initial point that satisfies all constraints
%and equality on the initial active set:

x0=[-3;4];

%Compute the gradient at x0: P*x0
g=P*x0;

%Positive-definite quadprog: (Eq. 5.47 in GMW)
pz=-(Z'*P*Z)\(Z'*g);
p0=Z*pz; %search direction

%Now we have to compute the maximum step length
%along the search direction (unity or less)
%Unity gives the exact step needed to reach the minimum
%subject to the equality constraint being considered

%Take the unity step and check constraints:
xtest=x0+p0;
A*xtest-B
%Indicates that constraints 1 and 3 would be violated.

%Find the maximum feasible step in that direction
alph_max1=(B(1)-A(1,:)*x0)/(A(1,:)*p0);
alph_max3=(B(3)-A(3,:)*x0)/(A(3,:)*p0);
alph_max=min(alph_max1,alph_max3);

%Find a step length not larger than alph_max that reduces
%the objective function. Here it's equal to alph_max

%Take that step
x1=x0+alph_max*p0;

%Because we had shortened steps for constraints 1 and 3,

%Find Lagrange multiplier *estimates* at x0 under current
%active set:

%A'*lambda=gtest

gtest=P*x1; 
lambda2=-3/4; %Negative: indicates this constraint can be 
%deleted (this may be incorrect, because it's based on an active
%set guess

%Add constraint 1, delete constraint 2,

A1=A(1,:);
%Find the null space basis
Z=null(A1);
%x1 is obviously still feasible

%Compute the gradient at x1: P*x1
g=P*x1;

%Positive-definite quadprog: (Eq. 5.47 in GMW)
pz=-(Z'*P*Z)\(Z'*g);
p1=Z*pz; %search direction

%Now we have to compute the maximum step length
%along the search direction (unity or less)
%Unity gives the exact step needed to reach the minimum
%subject to the equality constraint being considered

%Take the unity step and check constraints:
xtest=x1+p1;
A*xtest-B
%Indicates that no constraints would be violated.

%The maximum feasible step in that direction is alpha_max=1
alph_max=1;

%Find a step length not larger than alph_max that reduces
%the objective function. Here it's equal to alph_max

%Take that step
x2=x1+alph_max*p1;

%In quadratic programming, this step is the exact one to the
%minimum of the active set constrained problem

%We know we have reached the solution 

%For verification 
%Find Lagrange multipliers at x2 under current
%active set:

%A'*lambda=gtest
gtest=P*x2; 
lambda2=0.3636; %non-negative, so we have the solution

%Add some embellishments to the figure

%Shade the feasible set
%Vertices:[-4.5 5], [-2 0] [0 2])

patch([-4.5 -2 0],[5 0 2],'y')

%Show the iterations
plot(x0(1),x0(2),'*k')
plot(x1(1),x1(2),'*k')
plot(x2(1),x2(2),'*k')

title('Active set method example')
xlabel('x_1')
ylabel('x_2')







