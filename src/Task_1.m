close all;
clearvars;
tFormatStr = '%.3f';    %формат вывода времени на графике
V0 = 1;                 %продольная скорость движения пластины
c = sqrt(2);            %C - скорость распространения волны
dx = 0.01;              %h - шаг по длинне
dt = 0.01;              %tau - шаг по времени
Lx = 1;                 %длинна между роликами
D = 1;                  %? - неизвестный параметр
S = Lx;                 %? - длина пластины(абсолютная?)
Ro = 1;                 %плотность материала пластины
m = Ro*S;               %масса на единицу длины пластины
T0 = 1;                 %продольное напряжение пластины (на единицу длины)
T = 1;                  %время (верхняя граница)
X = 0:dx:Lx;            %сетка по длине пластины
Nx = size(X, 2);        %число узлов в сетке X
g_max = 1;              %max управляющей функции
g_max = g_max/m;
val_g = 0;              %начальное управление
g = @(t, val) 0;        %управляющая функция

NewFigure('Прямая задача');

W = 0.2*sin(pi*X);	    %начальное распределение точек листа
V = zeros(1,Nx);        %начальное распределение скорости листа
U = -pi*pi*0.2*sin(pi*X);

P = [ W; V; U];
P(:,1) = 0;             %граничное левое
P(:,Nx) = 0;            %граничное правое

Alpha = zeros(3,3,Nx);
Beta = zeros(3,1,Nx);

nu = dt/(dx^2);
A = [ 
        0 0 0; 
        0 -V0*dt/dx 0; 
        1 0 0 
    ];
B = [ 
        1 -dt 0; 
        0 1 (V0^2-c^2)*dt; 
        -2 0 -dx^2 
    ];
C = [ 
        0 0 0; 
        0 V0*dt/dx 0; 
        1 0 0 
    ];

hold all;
%axis([0 Lx 0 0.2]);
grid on;
plot(0,0, 'w');
p1 = plot(X,P(1,:));
set(p1,'LineWidth', 4);
legend(['t=' num2str(0,tFormatStr)],'W');
set (p1, 'Xdata', X, 'Ydata', P(1,:));
drawnow;  
for t=0:dt:T-dt
    F = [ P(1,:); P(2,:) + g(t+dt,val_g)*dt; zeros(1,Nx) ];
    for i=3:Nx
        Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
        Beta(:,:,i) = (A*Alpha(:,:,i-1)+B)\(F(:,i-1)-A*Beta(:,:,i-1));
    end
    for i=Nx-1:-1:2
        P(:,i) = Alpha(:,:,i+1)*P(:,i+1)+Beta(:,:,i+1);
    end
    set (p1, 'Xdata', X, 'Ydata', P(1,:));
    legend(['t=' num2str(t+dt,tFormatStr)],'W');
    drawnow;  
end

% обратная задача
f = NewFigure('Обратная задача');
figure(f);
%axis([0 Lx -2.5 1.5]);
% q1 = zeros(1,Nx);
% q2 = zeros(1,Nx);
% for i = 2:Nx-1
%     q1(i) = q1(i) + Ro*V0*(P(2,i+1)-P(2,i-1))/(2*dx);
%     q2(i) = -Ro*(P(2,i)+V0*(P(1,i+1)-P(1,i-1))/(2*dx));
% end
q1 = P(3,:) .* P(3,:);
q2 = P(1,:) .* P(1,:);
% for i = 2:Nx-1
%     q1(i) = q1(1) + Ro*V0*P(2,i);
%     q2(i) = -Ro*(P(2,i)+V0*(P(1,i+1)-P(1,i-1))/(2*dx));
% end
Q = [q1; q2];
Q(:,1) = 0;
Q(:,Nx) = 0;
dt = -dt;
nu = dt/(dx^2);
A = [ 
        0 nu*(c^2-V0^2); 
        0 -dt/dx*V0
    ];
B = [ 
        1 -2*nu*(c^2-V0^2); 
        dt 1 
    ];
C = [ 
        0 nu*(c^2-V0^2); 
        0 dt/dx*V0
    ];
grid on;
hold all;
plot(0,0, 'w');
p5 = plot(X,Q(1,:));
p6 = plot(X,Q(2,:));
legend(['t=' num2str(T,tFormatStr)],'q1','q2');
drawnow;
Alpha = zeros(2,2,Nx);
Beta = zeros(2,1,Nx);
for t=T:dt:0-dt
    F = [ Q(1,:); Q(2,:)];
    for i=2:Nx
        Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
        Beta(:,:,i) = (A*Alpha(:,:,i-1)+B)\(F(:,i-1)-A*Beta(:,:,i-1));
    end
    for i=Nx-1:-1:2
        Q(:,i) = Alpha(:,:,i+1)*Q(:,i+1)+Beta(:,:,i+1);
    end
    set (p5, 'Xdata', X, 'Ydata', Q(1,:));
    set (p6, 'Xdata', X, 'Ydata', Q(2,:));
    legend(['t=' num2str(t+dt,tFormatStr)],'q1','q2');
    drawnow;
    pause(1);
end

