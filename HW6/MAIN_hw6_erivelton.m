clc
clear all
close all

disp('This will take some time, because I reduced the step size for better animation results.');

[car1, car2, car3, car4, car1s, car2s, car3s, car4s, tsim, T, dt1, dt2, dt3, dt4] = ...
    serialCarIntesection;

disp('DONE');

animate(car1s, car2s, car3s, car4s, tsim)

figure;
ax1 = subplot(4,4,[1:8]);  bar(T'); ylabel('Centralized Solution')
ax2 = subplot(4,4,[9:16]); bar([zeros(size(dt1)); dt1;dt2;dt3;dt4]','stacked'); 
set(gca,'YDir','reverse'); 
ylabel('Serial Solution');
legend('car 1','car 2','car 3','car 4');

figure; hold on;
plot(car1s(1,:), 'LineWidth',2); 
plot(car2s(1,:), 'LineWidth',2); 
plot(car3s(2,:), 'LineWidth',2); 
plot(car4s(2,:), 'LineWidth',2); 
ylabel('Position'); xlabel('Iteration');
legend('car1','car2','car3','car4');

disp('running question 2');

ddMPCquestion2