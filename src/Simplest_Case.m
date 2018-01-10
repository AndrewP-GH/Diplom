clearvars; clc; close all;
f = figure('NumberTitle','off', 'Units','normalized', 'OuterPosition',[0 0 1 1]);
    pause(0.00001);
    frame_h = get(f,'JavaFrame');
    set(frame_h,'Maximized',1); 
figure(f);
V_0 = 1;
C_0 = sqrt(2);
dx = 0.01;
dt = 0.01;
g = @(t) 0;
g1 = 0;
g2 = 0;
Lx = 1;
T0=1;
Ro=1;
tFormatStr = '%.3f';
T = 1;
X = 0:dx:Lx;
Nx = size(X,2);
i1 = fix(Nx/3);
i2 = fix(Nx*2/3);

Q = [ 0.2*sin(pi*X); zeros(1,Nx) ];
Q(1,1) = 0;
Q(1,Nx) = 0;
A = [ 
        0 0; 
        dt/(dx^2)*(V_0^2-C_0^2) -dt/dx*V_0
    ];
B = [ 
        1 -dt; 
        -2*dt/(dx^2)*(V_0^2-C_0^2) 1
    ];
C = [ 
        0 0; 
        dt/(dx^2)*(V_0^2-C_0^2) dt/dx*V_0 
    ];
Alpha = zeros(2,2,Nx);
Beta = zeros(2,1,Nx);
hold all;
axis([0 Lx -0.4 0.4]);
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
    F = [ Q(1,:); Q(2,:)];
    Q(2,i1)=Q(2,i1)+g1;
    Q(2,i2)=Q(2,i2)+g2;
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


f2 = NewFigure('Обратная задача');
figure(f2);

% axis([0 Lx -25 50]);
q1 = zeros(1,Nx);
q2 = zeros(1,Nx);
for i = 2:Nx-1
    q1(i) = V_0*((T0+Ro*V_0)*(Q(1,i+1)-2*Q(1,i)+Q(1,i+1))/(dx^2)+Ro*(Q(2,i+1)-Q(2,i-1))/(2*dx));
    q2(i) = -Ro*(V_0*(Q(1,i+1)-Q(1,i-1))/(2*dx)+Q(2,i));
end
Q = [q1; q2];
Q(2,1) = 0;
Q(2,Nx) = 0;
Alpha = zeros(2,2,Nx);
Beta = zeros(2,1,Nx);
dt = -dt;
Q(1,1) = 0;
Q(1,Nx) = 0;
A = [ 
        0 dt/(dx^2)*(C_0^2-V_0^2); 
        0 -dt/dx*V_0
    ];
B = [ 
        1 -2*dt/(dx^2)*(V_0^2-C_0^2); 
        dt 1
    ];
C = [ 
        0 dt/(dx^2)*(C_0^2-V_0^2); 
        0 dt/dx*V_0 
    ];
hold all;
p5 = plot(X,q1);
p6 = plot(X,q2);
for t=T:dt:0-dt
    F = [ Q(1,:); Q(2,:)];
    for i=3:Nx
        Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
        Beta(:,1,i) = (A*Alpha(:,:,i-1)+B)\(F(:,i-1)-A*Beta(:,1,i-1));
    end
    for i=Nx-1:-1:2
        Q(:,i) = Alpha(:,:,i+1)*Q(:,i+1)+Beta(:,1,i+1);
    end
    set (p5, 'Xdata', X, 'Ydata', Q(1,:));
    set (p6, 'Xdata', X, 'Ydata', Q(2,:));
    legend(['t=' num2str(t+dt,tFormatStr)],'q1','q2');
    drawnow;
end