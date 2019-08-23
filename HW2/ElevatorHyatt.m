% Erivelton Gualter

clear all
close all

% Select Data 
% load('From_LL3-22_floor.mat');
% load('From-LL3-to-22.mat');
load('From-22-to-LL3.mat');

% Acceleration 
accX = Acceleration.X;
accY = Acceleration.Y;
accZ = Acceleration.Z;
acc = [accX accY accZ];

time = datevec(Acceleration.Timestamp);

initTime = time(1,5)*60 + time(1,6);
finalTime = time(end,5)*60 + time(end,6);

duration = finalTime - initTime;
Ts = duration/(length(time)-1);

tacel = 0:Ts:duration;

subplot(211); plot(tacel, acc(:,1), tacel, acc(:,2), tacel, acc(:,3)); 
    title('Acceleration Profile'); 
    ylabel('Accel X, Y, and Z'); 
    xlim([min(tacel) max(tacel)]);
    
subplot(212); plot(tacel, sqrt(acc(:,1).^2 + acc(:,2).^2 + acc(:,3).^2)); 
    ylabel('Absolute Accel'); 
    xlim([min(tacel) max(tacel)]);
    xlabel('Time [s]');
