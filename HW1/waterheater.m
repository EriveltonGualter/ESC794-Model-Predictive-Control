% Erivelton Gualter 
% 09/14/2018

% Dynamic Matrix Control Example

% Based on notes: Sergio Andres Castaño Giraldo
% https://controlautomaticoeducacion.com

% close all
clear all

% DMC parameters
P = 10; % Prediction Horizon
m = 5;  % Control Horizon
Q = 10*eye(P);  % States 
R = 1;          % Control

% Model
T=0.5;
z = tf('z');
Gp = 0.2731*z^(-3) / (1-0.8351*z^(-1));  
gi = step(Gp); 
Nm = length(gi)-1;
[Bp,Ap]=tfdata(Gp,'v');   % Get num and den
d = 1;

% Disturbance
Gd = 0.05*z^(-3) / (1-0.9*z^(-1));
gd = step(Gd); 
[Bd,Ad]=tfdata(Gd,'v');   % Get num and den

% Dynamic Matrix G
G = zeros(P,m); 
G(:,1)=gi(1:P);

for i=2:m
    for j=2:P
        G(j,i) = G(j-1,i-1);
    end
end

K = inv(G'*Q*G + R*eye*(Nm))*G';
K1 = K(1,:);

% Simulation

% Refereence
R = [ones(30,1);zeros(30,1);ones(30,1);zeros(30,1)];
% r = R;
% plot(R);

% Init.
nit = 120;
inc_u = 0;
u_ant(1:10) = 0;
u(1:10) = 0; ym(1:10) = 0; r(1:10) = 0;
% Referência
 
r(9:80) = 2; r(81:140) = 3; r(141:180) = 2; r(181:nit) = 2;
% r = R;
do(1:179) = 0;do(180:nit) = 0.1;
 
duf = zeros(1,Nm); %Acao de controle livre (Delta U Free)
 
w = 0; 

for k=m+1:nit-w
    ym(k) = Bp(1)*u(k-d)+Bp(2)*u(k-1-d)+Bp(3)*u(k-2-d)+Bp(4)*u(k-3-d)+Bp(5)*u(k-4-d)-Ap(2)*ym(k-1)-Ap(3)*ym(k-2)-Ap(4)*ym(k-3)+do(k);

    f = zeros(1,P);
    for kk=1:P
        for i=1:Nm-P
            vect_g(i)=gi(kk+i)-gi(i);
        end
        for i=Nm-P+1:Nm
            vect_g(i)=gi(Nm)-gi(i);
        end
        f(kk)=ym(k)+vect_g*duf';
    end
    
    inc_u=K1*(r(k+w)-f');

    if k==1
        u(k)=inc_u;
    else
        u(k)=u(k-1)+ inc_u;
    end

    aux_u=u_ant(1:length(Bp)-1);
    u_ant=[u(k) aux_u];

    aux_2=duf(1:Nm-1);
    duf=[inc_u aux_2];

end

nm=nit;
t = 0:T:(nm-1)*T;
figure
subplot(2,1,1); plot(t(1:nit-w),r(1:nit-w),'--k',t(1:nit-w),ym,'-r','Linewidth',2); xlabel('Time [s]'); ylabel('Responce'); legend('y_r','y','Location','SouthEast'); grid on; hold;
subplot(2,1,2); plot(t(1:nit-w),u,'b','Linewidth',2); xlabel('Time [s]'); ylabel('Control Input'); legend('u'); grid on;
