clearvars; clc; close all;
f = figure('NumberTitle','off', 'Units','normalized', 'OuterPosition',[0 0 1 1]);
    pause(0.00001);
    frame_h = get(f,'JavaFrame');
    set(frame_h,'Maximized',1); 
figure(f);
V_0 = 1;
C = sqrt(2);
dx = 0.01;
dt = 0.001;
g = @(t) 0;
Lx = 1;
T = 10;
X = 0:dx:Lx;
Nx = size(X,2);
Q = [ sin(pi*X); zeros(1,Nx) ];
Q(1,1) = 0;
Q(1,Nx) = 0;
A = [ 
        0 0; 
        dt/(dx^2)*(V_0^2-C^2) -dt/dx*V_0
    ];
B = [ 
        1 -dt; 
        -2*dt/(dx^2)*(V_0^2-C^2) 1
    ];
C = [ 
        0 0; 
        dt/(dx^2)*(V_0^2-C^2) dt/dx*V_0 
    ];
Alpha = zeros(2,2,Nx);
Beta = zeros(2,1,Nx);
hold all;
axis([0 Lx -1 1]);
p1 = plot(X,Q(1,:));
plot(0,0);plot(0,0);plot(0,0);plot(0,0);
Min = min(Q(1,:));
Max = max(Q(1,:));
localMin = Min;
localMax = Max;
legend(['t = ' num2str(0)]);
title({
    ['Min = ' num2str(Min)]
    ['localMin = ' num2str(localMin)]
    ['Max = ' num2str(Max)]
    ['localMax = ' num2str(localMax)]})
for t=0:dt:T-dt
    F = [ Q(1,:); Q(2,:)+dt*g(t) ];
    Alpha(:,:,3) = -B\C;
    Beta(:,:,3) = B\F(:,2);
    for i=4:Nx
        Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
        Beta(:,:,i) = (A*Alpha(:,:,i-1)+B)\(F(:,i-1)-A*Beta(:,:,i-1));
    end
    for i=Nx-1:-1:2
        Q(:,i) = Alpha(:,:,i+1)*Q(:,i+1)+Beta(:,:,i+1);
    end
    set (p1, 'Xdata', X, 'Ydata', Q(1,:));
    drawnow;
    legend(['t = ' num2str(t+dt)]);
    localMin = min(Q(1,:));
    localMax = max(Q(1,:));
    if localMin < Min
        Min = localMin;
    end
    if localMax > Max
        Max = localMax;
    end
    title({
        ['Min = ' num2str(Min)]
        ['localMin = ' num2str(localMin)]
        ['Max = ' num2str(Max)]
        ['localMax = ' num2str(localMax)]})
end