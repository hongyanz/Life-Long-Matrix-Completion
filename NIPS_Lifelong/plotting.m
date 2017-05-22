clc;
clear;

load('50times500passive.mat')
success_passive = success;
load('100times1000.mat')
success = success/10+success_passive;
success=success/2;
imshow(success);
axis on;
set(gca, 'xticklabel', [0.2:0.2:1]);
set(gca, 'yticklabel', [0.1:0.1:1]);
xlabel('Rank/m', 'fontsize', 20);
ylabel('Observations/m', 'fontsize', 20);
title('50\times500','fontsize', 20)