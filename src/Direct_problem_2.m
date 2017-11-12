clearvars;
tFormatStr = '%.3f';
V0 = 1;
c = sqrt(2);
T0 = 1;
dx = 0.01;
dt = 0.01;
Lx = 1;
D = 1;
S = Lx;
ro = 1;                             %Ro
m = ro*S;
val_g = 0;
g = @(t,val) 0;
g_max = 1/m;
T = 1;
X = 0:dx:Lx;
Nx = size(X,2);
littel_beta = 0.1;

filename = 'testnew51.gif';
gifDelayTime=0;

f = NewFigure('Прямая задача');
figure(f);
W = 0.2*sin(pi*X); %zeros(1,Nx);    %начальное распределение точек листа
V = zeros(1,Nx); %zeros(1,Nx);      %начальное распределение скорости листа
U = -pi^2*0.2*sin(pi*X);
P = [ W; V; U];
P(:,1) = 0;                         %граничное левое
P(:,Nx) = 0;                        %граничное правое

Alpha = zeros(3,3,Nx);
Beta = zeros(3,1,Nx);

nu = dt/(dx^2);
A = [ 
        0 0 0; 
        nu*(V0^2-c^2) -V0*dt/dx D/m*nu; 
        1 0 0 
    ];
B = [ 
        1 -dt 0; 
        -2*nu*(V0^2-c^2) 1 -2*D/m*nu; 
        -2 0 -dx^2 
    ];
C = [ 
        0 0 0; 
        nu*(V0^2-c^2) V0*dt/dx D/m*nu; 
        1 0 0 
    ];

hold all;
%axis([0 Lx 0 0.2]);
grid on;
p0 = plot(0,0, 'w');
p1 = plot(X,P(1,:));
set(p1,'LineWidth', 4);
legend(['t=' num2str(0,tFormatStr)],'W');
set (p1, 'Xdata', X, 'Ydata', P(1,:));
drawnow;  
SaveAsGif(filename, gifDelayTime, 0);
for t=0:dt:T-dt
    F = [ P(1,:); P(2,:) + g(t,val_g)*dt; zeros(1,Nx) ];
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
    SaveAsGif(filename, gifDelayTime, 1);
end

% обратная задача
f = NewFigure('Обратная задача');
figure(f);
%axis([0 Lx -1 1]);
q1 = P(3,:)*(T0+ro*V0^2)+ ro*V0*P(2,:);
q2 = zeros(1,Nx);
q2(1) = -ro*(P(2,1)+V0*(P(1,2)-P(1,1))/dx);
q2(Nx) = -ro*(P(2,Nx)+V0*(P(1,Nx)-P(1,Nx-1))/dx);
for i = 2:Nx-1
    q2(i) = -ro*(P(2,i)+V0*(P(1,i+1)-P(1,i-1))/(2*dx));
end
Q = [q1; q2; zeros(1,Nx)];
Q(:,1) = 0;
Q(:,Nx) = 0;
Alpha = zeros(3,3,Nx);
Beta = zeros(3,1,Nx);
dt = -dt;
nu = dt/(dx^2);
A = [ 
        0 nu*(c^2-V0^2) -nu*D/m; 
        0 -nu*dx*V0 0;                      
        0 1 0 
    ];
B = [ 
        1 -2*nu*(c^2-V0^2) 2*nu*D/m; 
        dt 1 0; 
        0 -2 -dx^2 
    ];
C = [ 
        0 nu*(c^2-V0^2) -nu*D/m; 
        0 nu*dx*V0 0;                       
        0 1 0 
    ];
grid on;
hold all;
p4 = plot(0,0, 'w');
p5 = plot(X,Q(1,:));
p6 = plot(X,Q(2,:));
%set (p1,'LineWidth', 4);
legend(['t=' num2str(T,tFormatStr)],'q1','q2');
drawnow;
pause(5);
for t=T:dt:0-dt
    F = [ Q(1,:); Q(2,:); zeros(1,Nx) ];
    for i=3:Nx
        Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
        Beta(:,:,i) = (A*Alpha(:,:,i-1)+B)\(F(:,i-1)-A*Beta(:,:,i-1));
    end
    %Alpha(:,:,Nx) = zeros(3,3);
    %Beta(:,:,Nx) = zeros(3,1);
    for i=Nx-1:-1:2
        Q(:,i) = Alpha(:,:,i+1)*Q(:,i+1)+Beta(:,:,i+1);
    end
    set (p5, 'Xdata', X, 'Ydata', Q(1,:));
    set (p6, 'Xdata', X, 'Ydata', Q(2,:));
    legend(['t=' num2str(t+dt,tFormatStr)],'q1','q2');
%     pause(1);
    drawnow;
end

temp_g=0;
for i=1:1:Nx-1
    temp_g = temp_g + (Q(2,i)+Q(2,i+1))/2*dx;
end
temp_g=temp_g/(2*littel_beta);
if (temp_g<0)
    val_g = 0;
else
    if (temp_g>g_max)
        val_g = g_max;
    else
        val_g = temp_g;
    end
end