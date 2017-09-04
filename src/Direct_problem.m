clearvars; clc;
tFormatStr = '%.3f';
V_0 = 1;
C = sqrt(2);
dx = 0.05;
dt = 0.01;
m = 1;
g = @(t)t/m;
Lx = 1;
T = 10;
X = 0:dx:Lx;
Nx = size(X,2);
Q = [ 0.2*sin(pi*X); zeros(1,Nx) ];
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
p1 = plot(X,Q(1,:));
%p2 = plot(X,Q(2,:));
legend(['t=' num2str(0,tFormatStr)],'W','V');
for t=0:dt:T-dt
    D = [ Q(1,:); Q(2,:)+dt*g(t) ];
    for i=3:Nx
        Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
        Beta(:,:,i) = (A*Alpha(:,:,i-1)+B)\(D(:,i-1)-A*Beta(:,:,i-1));
    end
    for i=Nx-1:-1:2
        Q(:,i) = Alpha(:,:,i+1)*Q(:,i+1)+Beta(:,:,i+1);
    end
    set (p1, 'Xdata', X, 'Ydata', Q(1,:));
   % set (p2, 'Xdata', X, 'Ydata', Q(2,:));
    drawnow;
    legend(['t=' num2str(t+dt,tFormatStr)],'W','V');
end