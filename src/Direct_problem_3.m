clearvars;
tFormatStr = '%.3f';
V0 = 1;
c = sqrt(2);
T0 = 1;
dx = 0.005;
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

f2 = NewFigure('Прямая задача');
figure(f2);
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
legend(['t=' num2str(0,tFormatStr)]);

q1 = P(3,:)*(T0+ro*V0^2)+ ro*V0*P(2,:);
    %q2 = zeros(1,Nx);
q2(1) = -ro*(P(2,1)+V0*(P(1,2)-P(1,1))/dx);
q2(Nx) = -ro*(P(2,Nx)+V0*(P(1,Nx)-P(1,Nx-1))/dx);
for i = 2:Nx-1
    q2(i) = -ro*(P(2,i)+V0*(P(1,i+1)-P(1,i-1))/(2*dx));
end
%p1 = plot(X,q1);
p2 = plot(X,q2);


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
    
    q1 = P(3,:)*(T0+ro*V0^2)+ ro*V0*P(2,:);
    %q2 = zeros(1,Nx);
    q2(1) = -ro*(P(2,1)+V0*(P(1,2)-P(1,1))/dx);
    q2(Nx) = -ro*(P(2,Nx)+V0*(P(1,Nx)-P(1,Nx-1))/dx);
    for i = 2:Nx-1
        q2(i) = -ro*(P(2,i)+V0*(P(1,i+1)-P(1,i-1))/(2*dx));
    end

    %set (p1, 'Xdata', X, 'Ydata', q1);
    set (p2, 'Xdata', X, 'Ydata', q2);
    legend(['t=' num2str(t+dt,tFormatStr)],'W');
    
    drawnow;  
    SaveAsGif(filename, gifDelayTime, 1);
end