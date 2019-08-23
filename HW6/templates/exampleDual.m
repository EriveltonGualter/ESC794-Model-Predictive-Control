%exampleDual

%Progressive example of solving a quadratic program by duality methods and
%decomposition

%Programmed by H Richter for MPC course: Selected Topics on Engineering
%Science

%Cleveland State University, Mechanical Engineering Department
%Fall 2018


%Minimize F(x)=0.5*x'*P*x+c'*x

%subject to Ax<=B

%Some Pdef matrix:

P=diag([1 2 4]);    

A=-[1 1 1; ...
    0 0 -1;...
    0 -1 0;...
    -1 0 0];

B=-[4;-2;0;-4];
c=[1;-1;0];

%Find quadprog solution: primal problem
[Xprimal,Fval]=quadprog(P,c,A,B)


%Formulate the dual problem: 
%min 0.5*lambda'*A*inv(P)*A*lambda+c'*inv(P)*A'*lambda
%subject to lambda_i >=0

H=A*(P\A');
f=A*(P\c)+B;

Adual=-eye(4);Bdual=[0;0;0;0];

[Lambda,Dval]=quadprog(H,f,Adual,Bdual)

%Correct Dval to add terms that were indepedent of Lambda
Dval=-(Dval+c'*(P\c))

%Check that Lambda can be used to find the primal solution
XprimalCheck=-P\(A'*Lambda+c)


%Now use an iterative scheme, first with a decoupled problem (remove first
%inequality constraint

%The primal inequalities are

% x1+x2<= 2
% x3   <= 4

%Partitioning: Note the cost function is already partitioned
%{x1 x2},{x3}

%Matrices
A=[1 1 0; ...
   0 0 1]; 
B=[2;4];

%Find quadprog solution: primal problem
[Xprimal,Fval]=quadprog(P,c,[],[],A,B)

%Enter distributed optimization
%Initialize Lagrange multiplier vector 
Lambda=[0;0];

%Solve primal problems *in parallel* for the subystems

A12=A(:,1:2);
A3=A(:,3);

P12=diag([1 2]); P3=4;
c12=[1;-1];c3=0;
rho=0.99; %convergence factor


for i=1:25
%{1,2}: %Minimize Xprimal with current Lambda
X12=-P12\(A12'*Lambda+c12);
%{3}
X3=-P3\(A3'*Lambda+c3);

%Reconstruct primal solution
X=[X12;X3];

%Update multiplier
Lambda=Lambda+rho*(A*X-B);
end
    