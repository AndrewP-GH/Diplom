clearvars; close all;
addpath('./Functions');
Figure = NewFigure('������');

%% ����� ��������� ������
global FormatStr        %������ ������ ������� �� �������
FormatStr = '%.3f';
V_0 = 0.5;              %���������� �������� �������� ��������
C = sqrt(2);            %C - �������� ��������������� �����
V_C = V_0^2 - C^2;      
global dx;
global dt;
dx = 0.01;              %h - ��� �� ������
dt = 0.01;              %tau - ��� �� �������
s_dt = Significant(dt);
L_x = 1;                %������ ����� ��������
S = L_x;                %? - ����� ��������(����������?)
Ro = 1;                 %��������� ��������� ��������
m = Ro*S;               %����� �� ������� ����� ��������
T_0 = 1;                %���������� ���������� �������� (�� ������� �����)

T = Roundx(2*L_x*C/(-V_C), s_dt, 'ceil');                  %����� (������� �������)
Tarray = 0:dt:T;
Nt = length(Tarray);    %����� ����� ���������

X = 0:dx:L_x;           %����� �� ����� ��������
Nx = size(X, 2);        %����� ����� � ����� X
g_max = 1;              %max ����������� �������
i_1 = round(Nx/3);
i_2 = round(Nx/3*2);
iterations = 2;         %����� ��������
Ymax = 1;

global CalcExtrems      %�������� ������� � ��������
CalcExtrems = true;
OnlyQ2 = true;

folder = CreateImageFolder([datestr(now, 'dd-mmm-yyyy HH_MM_SS') '_������']);
image_type = '.tiff';
gif_type = '.gif';
image_name = 'tmp';
save_gif = false;
gif_delay = 1/24;

%% ��������� ������ ������
P_0 = 0.2*sin(pi*X);	%��������� ������������� ����� �����
P_1 = zeros(1,Nx);      %��������� ������������� �������� �����
nu = dt / dx^2;
mu = dt / dx;
A_i =   [   0           0; 
            nu*V_C      -mu*V_0
        ];
B_i =   [   1           -dt;
            -2*nu*V_C   1
        ];
C_i =   [   0           0; 
            nu*V_C      mu*V_0
        ];
L = zeros(2,2,Nx);      %����������� ������������ ������ ������
M = zeros(2,1,Nx);      %����������� ������������ ������ ������

L(:,:,3) = - B_i \ C_i;
for i = 4:Nx
    L(:,:,i) = -(A_i * L(:,:,i-1) + B_i) \ C_i;
end

%% ��������� �������� ������
dt_ = -dt;
nu_ = dt_ / dx^2;
mu_ = dt_ / dx;
A_i_ =  [   0	-nu_*V_C;
            0	-mu_*V_0
        ];
B_i_ =  [   1	-2*nu_*(-V_C);
            dt_	1
        ];
C_i_ =  [   0	-nu_*V_C;
            0	mu_*V_0
        ];
L_ = zeros(2,2,Nx);     %����������� ������������ �������� ������
M_ = zeros(2,1,Nx);     %����������� ������������ �������� ������
Q = zeros(2,Nx);        %������ �������� (q_1, q_2)

L_(:,:,2) = [   0   (2*nu_*V_0*V_C) / (V_0 - mu_*V_C);  %��������� ������� � �������� ������� (dt_<0), �� ��������� � V
                0   0
            ];
for i = 3:Nx
    L_(:,:,i) = -(A_i_ * L_(:,:,i-1) + B_i_) \ C_i_;
end
Ew = zeros(1,iterations);
QH = zeros(Nt,Nx);      %��� ����� ��������� �������� q_2 � ������ ������ �������
Control = zeros(Nt,Nx);

%% ���� ����������
for k=1:iterations
    %% ���������� ������ ������
    [p_1, p_2] = FigurePrepare(Figure, X);
    axis([0 L_x -Ymax Ymax]);
    p_1.LineWidth = 3;
    P = [P_0; P_1];         %������ �������� (W, V)
    P(:,1) = 0;             %��������� �����
    P(:,Nx) = 0;            %��������� ������
    SetTwoLinesInPlots(p_1, p_2, P, 0, 'W', 'V');
    if CalcExtrems == true
        [Min, Max] = LocalExtrems(P, 0, 0, 1, 'W');
    end
    image_name = num2str(k);
    if save_gif == true
        SaveAsGif(folder, [image_name gif_type], gif_delay, 0);
    end
    SaveAsGif(folder, [image_name '_begin' image_type], 1, 0);
    for t = 1:Nt-1
        F = [   P(1,:); 
                P(2,:)
            ];
        x1 = 1 : round(Nx/2);
        x2 = round(Nx/2)+1 : Nx;
        integral = trapz(X(x1), QH(t,x1));
        if integral > 0
            Control(t, x1) = g_max;
        end
        integral = trapz(X(x2), QH(t,x2));
        if integral < 0
            Control(t, x2) = -g_max;
        end 
        F(2,:) = F(2,:) + dt*Control(t,:);
              
        M(:,:,3) = B_i \ F(:,2);
        for i = 4:Nx
            M(:,:,i) = -(A_i * L(:,:,i-1) + B_i) \ (A_i * M(:,:,i-1) - F(:,i-1));
        end
        for i = Nx-1:-1:2
            P(:,i) = L(:,:,i+1) * P(:,i+1) + M(:,:,i+1);
        end
        SetTwoLinesInPlots(p_1, p_2, P, t*dt, 'W', 'V');
        if CalcExtrems == true
            [Min, Max] = LocalExtrems(P, Min, Max, 1, 'W');
        end
        if save_gif == true
            SaveAsGif(folder, [image_name gif_type], gif_delay, 1);
        end
    end
    Ew(k) = PaperFullEnergy(P, T_0, Ro, V_0);
    disp([ 'Full energy Ew = ' num2str(Ew(k)) ]);
    SaveAsGif(folder, [image_name '_end' image_type], 1, 0);
    %% ���������� �������� ������
    if k ~= iterations
        [p_1, p_2] = FigurePrepare(Figure, X);
        axis([0 L_x -Ymax Ymax]);
        W_xx = (P(1,1:Nx-2) - 2*P(1,2:Nx-1) + P(1,3:Nx))/dx^2;
        Q(1,2:Nx-1) = T_0*W_xx ...
            + Ro*V_0 * ( ...
                (P(2,3:Nx) - P(2,1:Nx-2))/(2*dx) ...
                + V_0*W_xx ...                                  
            );
        Q(2,2:Nx-1) = -Ro * ( ...
                P(2,2:Nx-1) + V_0*(P(1,3:Nx) - P(1,1:Nx-2))/(2*dx) ...                              
            );
        W_xx = (2*P(1,1) - 5*P(1,2) + 4*P(1,3) - P(1,4))/dx^2;
        Q(1:1) = T_0*W_xx ...
            + Ro*V_0 * ( ...
                (-3/2*P(2,1) + 2*P(2,2) - 1/2*P(2,3))/dx ...
                + V_0*W_xx ...                                  
            );
        Q(2,1) = 0;
        W_xx = (2*P(1,Nx) - 5*P(1,Nx-1) + 4*P(1,Nx-2) - P(1,Nx-3))/dx^2;
        Q(1,Nx) = T_0*W_xx ...
            + Ro*V_0 * ( ...
                (3/2*P(2,Nx) - 2*P(2,Nx-1) + 1/2*P(2,Nx-2))/dx ...
                + V_0*W_xx ...                                  
            );
        Q(2,Nx) = 0;
        if OnlyQ2 == true
            SetTwoLinesInPlots(p_1, p_2, Q, T, 'q2');
        else
            SetTwoLinesInPlots(p_1, p_2, Q, T, 'q1', 'q2');
        end
        if CalcExtrems == true
            [Min, Max] = LocalExtrems(Q, 0, 0, 2, 'q2');
        end
        QH(Nt,:) = Q(2,:);
        image_name = ['q_' num2str(k)];
        if save_gif == true
            SaveAsGif(folder, [image_name gif_type], gif_delay, 0);
        end
        for t=Nt-1:-1:1
            F = [   Q(1,:); 
                    Q(2,:)
                ];
            M_(:,:,2) = [   (V_0 * Q(1,1)) / (V_0 - mu_*V_C);
                            0
                        ];
            for i = 3:Nx
                M_(:,:,i) = -(A_i_ * L_(:,:,i-1) + B_i_) \ (A_i_ * M_(:,:,i-1) - F(:,i-1));
            end
            Q(1,Nx) = V_0 * (Q(1,Nx) + 2*M_(2,1,Nx)*nu_*V_C) / (V_0 + mu_*V_C(1 - 2*L_(2,1,Nx)*V_0)/dx);
            Q(2,Nx) = 0;
            for i = Nx-1:-1:1
                Q(:,i) = L_(:,:,i+1) * Q(:,i+1) + M_(:,:,i+1);
            end
            Q(2,1) = 0;
            if OnlyQ2 == true
            	SetTwoLinesInPlots(p_1, p_2, Q, (t-1)*dt, 'q2');
            else
                SetTwoLinesInPlots(p_1, p_2, Q, (t-1)*dt, 'q1', 'q2');
            end
            if CalcExtrems == true
                [Min, Max] = LocalExtrems(Q, Min, Max, 2, 'q2');
            end
            QH(t,:) = Q(2,:);
            if save_gif == true
                SaveAsGif(folder, [image_name gif_type], gif_delay, 1);
            end
        end
        SaveAsGif(folder, [image_name '_end' image_type], 1, 0);
    else
        FigurePrepare(Figure, Tarray, false);
        h = surf(X,Tarray, Control);
        set(h, 'LineWidth', 0.05)
        axis tight;
        box on;
%         shading interp
        xlabel('X');
        ylabel('T');
        zlabel('Control')      
        newmap = jet();  
        ncol = size(newmap, 1);
        newmap(2:ncol-1,:) = ones(ncol-2,3);
        newmap(1,:) = [0 0 1];
        newmap(ncol,:) = [1 0 0];
        colormap(newmap);
        colorbar;
        view(0,90);
        title({
            strcat('Ew_{���.} =', [' ' num2str(Ew(1))])
            strcat('Ew_{���.} =', [' ' num2str(Ew(k))])
        });
        SaveAsGif(folder, ['control' image_type], 1, 0);
    end
end