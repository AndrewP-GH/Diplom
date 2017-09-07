clearvars; clc;
tFormatStr = '%.3f';
V_0 = 1;
C = sqrt(2);
dx = 0.001;
dt = 0.001;
D = 1;
m = 1;
g = @(t)0/m;
Lx = 1;
T = 1;
X = 0:dx:Lx;
Nx = size(X,2);

W = 0.2*sin(pi*X); %zeros(1,Nx);
V = zeros(1,Nx); %zeros(1,Nx);
U = -pi^2*0.2*sin(pi*X);
P = [ W; V; U];
P(:,Nx) = 0; %граничное

Alpha = zeros(3,3,Nx);
Beta = zeros(3,1,Nx);

nu = dt/(dx^2);
A = [ 
        0 0 0; 
        nu*(V_0^2-C^2) -V_0*dt/dx D/m*nu; 
        1 0 0 
    ];
B = [ 
        1 -dt 0; 
        -2*nu*(V_0^2-C^2) 1 -2*D/m*nu; 
        -2 0 -dx^2 
    ];
C = [ 
        0 0 0; 
        nu*(V_0^2-C^2) V_0*dt/dx D/m*nu; 
        1 0 0 
    ];

hold all;
axis([0 Lx -2.3 2.3]);
grid on;
p0 = plot(0,0, 'w');
p1 = plot(X,P(1,:));
p2 = plot(X,P(2,:));
%p3 = plot(X,P(3,:));
set(p1,'LineWidth', 4);
set(p2,'LineWidth', 2);
%set(p3,'LineWidth', 1);
legend(['t=' num2str(0,tFormatStr)],'W','V');
for t=0:dt:T-dt
    F = [ P(1,:); P(2,:) + g(t)*dt; zeros(1,Nx) ];
    for i=3:Nx
        Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
        Beta(:,:,i) = (A*Alpha(:,:,i-1)+B)\(F(:,i-1)-A*Beta(:,:,i-1));
    end
    for i=Nx-1:-1:2
        P(:,i) = Alpha(:,:,i+1)*P(:,i+1)+Beta(:,:,i+1);
    end
    set (p1, 'Xdata', X, 'Ydata', P(1,:));
    set (p2, 'Xdata', X, 'Ydata', P(2,:));
    %set (p3, 'Xdata', X, 'Ydata', P(3,:));
    legend(['t=' num2str(t+dt,tFormatStr)],'W','V');
    drawnow;
end