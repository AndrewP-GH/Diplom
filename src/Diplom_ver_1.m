clearvars; clc; close all;
f = figure('NumberTitle','off', 'Units','normalized', 'OuterPosition',[0 0 1 1]);
    pause(0.00001);
    frame_h = get(f,'JavaFrame');
    set(frame_h,'Maximized',1); 
figure(f);
Lx = 1;
V_0 = 1;
C = sqrt(2);
dx = 0.01;
dt = 0.001;
g = @(t) 0;