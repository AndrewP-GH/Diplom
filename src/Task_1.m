close all;
clearvars;
tFormatStr = '%.3f';    %������ ������ ������� �� �������
V0 = 1;                 %���������� �������� �������� ��������
c = sqrt(2);            %C - �������� ��������������� �����
dx = 0.01;              %h - ��� �� ������
dt = 0.01;              %tau - ��� �� �������
Lx = 1;                 %������ ����� ��������
D = 1;                  %? - ����������� ��������
S = Lx;                 %? - ����� ��������(����������?)
Ro = 1;                 %��������� ��������� ��������
m = Ro*S;               %����� �� ������� ����� ��������
T0 = 1;                 %���������� ���������� �������� (�� ������� �����)
T = 1;                  %����� (������� �������)
X = 0:dx:Lx;            %����� �� ����� ��������
Nx = size(X, 2);        %����� ����� � ����� X
g_max = 1;              %max ����������� �������
g_max = g_max/m;
val_g = 0;              %��������� ����������
g = @(t, val) 0;        %����������� �������

f = NewFigure('������ ������');
saver =  CSaveAsGif();

W = 0.2*sin(pi*X);	    %��������� ������������� ����� �����
V = zeros(1,Nx);        %��������� ������������� �������� �����
U = -pi^2*0.2*sin(pi*X);

P = [ W; V; U];
P(:,1) = 0;             %��������� �����
P(:,Nx) = 0;            %��������� ������

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
axis([0 Lx 0 0.2]);
grid on;
plot(0,0, 'w');
p1 = plot(X,P(1,:));
set(p1,'LineWidth', 4);
legend(['t=' num2str(0,tFormatStr)],'W');
set (p1, 'Xdata', X, 'Ydata', P(1,:));
drawnow;  
for t=0:dt:T-dt
    F = [ P(1,:); P(2,:) + g(t+dt,val_g)*dt; zeros(1,Nx) ];
    Alpha(:,:,2) = [ 0 0 0; 0 0 0; 0 2/(dx^2) 0];
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

% �������� ������
f = NewFigure('�������� ������');
figure(f);
%axis([0 Lx -2.5 1.5]);
q1 = P(3,:)*(T0+Ro*V0^2);
q2 = zeros(1,Nx);
% for i = 2:Nx-1
%     q1(i) = q1(i) + Ro*V0*(P(2,i+1)-P(2,i-1))/(2*dx);
%     q2(i) = -Ro*(P(2,i)+V0*(P(1,i+1)-P(1,i-1))/(2*dx));
% end
for i = 2:Nx-1
    q1(i) = q1(1) + Ro*V0*P(2,i)/dx;
    q2(i) = -Ro*(P(2,i)+V0*(P(1,i+1)-P(1,i-1))/(2*dx));
end
Q = [q1; q2; zeros(1,Nx)];
Q(:,1) = 0;
Q(:,Nx) = 0;
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
plot(0,0, 'w');
p5 = plot(X,Q(1,:));
p6 = plot(X,Q(2,:));
legend(['t=' num2str(T,tFormatStr)],'q1','q2');
drawnow;
pause(5);
tmpCoef = [ 
        0 0 -nu*D/m; 
        0 0 0;                       
        0 0 0 
    ];
Alpha = zeros(3,3,Nx);
Beta = zeros(3,1,Nx);
m = ceil(Nx/2);
for t=T:dt:0-dt
    F = [ Q(1,:); Q(2,:); zeros(1,Nx) ];
%     Alpha(:,:,2) = -tmpCoef;
%     for i=3:m-1
%         Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
%         Beta(:,:,i) = (A*Alpha(:,:,i-1)+B)\(F(:,i-1)-A*Beta(:,:,i-1));
%     end
%     Alpha(:,:,Nx-1) = -tmpCoef\A;
%     for i=Nx-2:-1:m
%         Alpha(:,:,i) = -(A*Alpha(:,:,i+1)+B)\C;
%         Beta(:,:,i) = (A*Alpha(:,:,i+1)+B)\(F(:,i)-A*Beta(:,:,i+1));
%     end
%     Q(:,m-1) = (1-Alpha(:,:,m)*Alpha(:,:,m))\(Alpha(:,:,m)*Beta(:,:,m)+Beta(:,:,m));
%     Q(:,m) = (1-Alpha(:,:,m)*Alpha(:,:,m))\(Alpha(:,:,m)*Beta(:,:,m)+Beta(:,:,m))*Alpha(:,:,m)+Beta(:,:,m);
%     for i=m-2:-1:2
%         Q(:,i) = Alpha(:,:,i+1)*Q(:,i+1)+Beta(:,:,i+1);
%     end
%     for i=m+1:Nx-1
%         Q(:,i) = Alpha(:,:,i)*Q(:,i-1)+Beta(:,:,i);
%     end
    for i=3:Nx
        Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
        Beta(:,:,i) = (A*Alpha(:,:,i-1)+B)\(F(:,i-1)-A*Beta(:,:,i-1));
    end
    for i=Nx-1:-1:2
        Q(:,i) = Alpha(:,:,i+1)*Q(:,i+1)+Beta(:,:,i+1);
    end
    set (p5, 'Xdata', X, 'Ydata', Q(1,:));
    set (p6, 'Xdata', X, 'Ydata', Q(2,:));
    legend(['t=' num2str(t+dt,tFormatStr)],'q1','q2');
    pause(1);
    drawnow;
end

