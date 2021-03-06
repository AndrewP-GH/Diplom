clearvars; close all;
addpath('./Functions');
Figure = NewFigure('������');
warning off;

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
L_x = 1;                %������ ����� ��������
D = 1;                  %? - ����������� ��������
S = L_x;                %? - ����� ��������(����������?)
Ro = 1;                 %��������� ��������� ��������
m = Ro*S;               %����� �� ������� ����� ��������
T_0 = 1;                %���������� ���������� �������� (�� ������� �����)

T = 1;                  %����� (������� �������)
Tarray = 0:dt:T;
Nt = length(Tarray);    %����� ����� ���������

X = 0:dx:L_x;           %����� �� ����� ��������
Nx = size(X, 2);        %����� ����� � ����� X
g = 0;
g_min = 0;
g_max = 200;              %max ����������� �������
h_g = 10;
Ew_max = 0.0001;
Eps_Ew = 0.0001;
i_1 = round(Nx/3);
i_2 = round(Nx/3*2);
iterations = 2;         %����� ��������
Ymax = 2;

global CalcExtrems      %�������� ������� � ��������
CalcExtrems = true;
OnlyQ2 = true;

folder = CreateImageFolder([datestr(now, 'dd-mmm-yyyy HH_MM_SS') '_������']);
image_type = '.tiff';
gif_type = '.gif';
image_name = 'tmp';
save_gif = true;
gif_delay = 1/24;

%% ��������� ������ ������
W = 0.2*sin(pi*X);	    %��������� ������������� ����� �����
V = zeros(1,Nx);        %��������� ������������� �������� �����
U = -pi*pi*0.2*sin(pi*X);
nu = dt / dx^2;
mu = dt / dx;
A_i =   [   0       0       0;
            nu*V_C  -mu*V_0 nu*D/m;
            1       0       0
        ];
B_i =   [   1           -dt	0;
            -2*nu*V_C   1	-2*nu*D/m;
            -2          0   -dx^2
        ];
C_i =   [   0       0       0;
            nu*V_C  mu*V_0  nu*D/m;
            1       0       0
        ];
L = zeros(3,3,Nx);      %����������� ������������ ������ ������
M = zeros(3,1,Nx);      %����������� ������������ ������ ������

L(:,:,3) = - B_i \ C_i;
for i = 4:Nx
    L(:,:,i) = -(A_i * L(:,:,i-1) + B_i) \ C_i;
end


%% ��������� �������� ������
dt_ = -dt;
nu_ = dt_ / dx^2;
mu_ = dt_ / dx;
A_i_ =  [   0	nu_*(-V_C)	-nu_*D/m;
            0	-mu_*V_0	0;
            0	1           0
        ];
B_i_ =  [   1   2*nu_*V_C	2*nu_*D/m;
            dt_	1           0;
            0   -2          -dx^2
        ];
C_i_ =  [   0	nu_*(-V_C)	-nu_*D/m;
            0	mu_*V_0     0;
            0	1           0
        ];
L_ = zeros(3,3,Nx);     %����������� ������������ �������� ������
M_ = zeros(3,1,Nx);     %����������� ������������ �������� ������
Q = zeros(3,Nx);        %������ �������� (q_1, q_2, q_3)

L_(:,:,2) = [   0 0      0;
                0 0      0;
                0 2/dx^2 0
            ];
for i = 3:Nx
    L_(:,:,i) = -(A_i_ * L_(:,:,i-1) + B_i_) \ C_i_;
end
Ew = zeros(1,iterations);
g_values = zeros(1,iterations);
QH = zeros(Nt,Nx);      %��� ����� ��������� �������� q_2 � ������ ������ �������
Control = zeros(Nt,Nx);
prevControl = zeros(Nt,Nx);

%% ���� ����������
k = 1;
iter_in_g = 0;      %����� �������� ������ ������ [g_min,g_max]
l2r = true;
while k ~= iterations+1
    %% ���������� ������ ������
    prevControl = Control;
    [p_1, p_2] = FigurePrepare(Figure, X);
    axis([0 L_x -Ymax Ymax]);
    p_1.LineWidth = 3;
    P = [W; V; U];         %������ �������� (W, V)
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
    for t = 1:Nt-1
        F = [   P(1,:);
                P(2,:);
                zeros(1,Nx)
            ];
        intervals = [ {1 : round(Nx/2)}; {round(Nx/2)+1 : Nx}];
        Control(t,:) = AddControlOnInterval(intervals, X, QH(t,:), -g);
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
    Ew(k) = PaperFullEnergy(P, T_0, Ro, V_0, X);
    g_values(k) = g;
    disp([ 'Full energy Ew = ' num2str(Ew(k)) '; g = ' num2str(g)]);
    iter_in_g = iter_in_g+1;
    if Ew(k) >= Ew_max && (k == 1 || abs( Ew(k) - Ew(k-1) ) >= Eps_Ew )% && h_g >= 0.001 && g ~= g_max && (g ~= g_min || k == 1)     %���� ��� ������� ��� ������� ������, � �� g==g_max ��� ������� ����������
        if (iter_in_g == 2 && Ew(k-1) < Ew(k)) || ...
                (iter_in_g > 2 && Ew(k-2) > Ew(k-1) && Ew(k-1) < Ew(k))     %��� ��� �/� g_min � g_max => ����� �������
            if l2r
                g_min = g - h_g;
                g_max = g;
                if iter_in_g > 2
                    g_min = g_min - h_g;
                end
            else
                g_min = g;
                g_max = g + h_g;
                if iter_in_g > 2
                    g_max = g_max + h_g;
                end
            end
            h_g = h_g / 5;
            l2r = ~l2r;
            iter_in_g = 0;
        end
        if h_g >= 0.01 && (iter_in_g == 0 || (g < g_max && g > g_min))
            iterations = iterations + 1;
        end
       	if l2r
            g = g + h_g;
        else
            g = g - h_g;
        end
    end
    SaveAsGif(folder, [image_name '_end' image_type], 1, 0);
    %% ���������� �������� ������
    if k == 1 %if k ~= iterations
        [p_1, p_2] = FigurePrepare(Figure, X);
        axis([0 L_x -Ymax Ymax]);
        W_xx = (P(1,1:Nx-2) - 2*P(1,2:Nx-1) + P(1,3:Nx))/dx^2;
        Q(1,2:Nx-1) = (T_0 + Ro*V_0^2) * W_xx ...
            + Ro*V_0 * (P(2,3:Nx) - P(2,1:Nx-2))/(2*dx);
        Q(2,2:Nx-1) = -Ro * ( ...
                P(2,2:Nx-1) + V_0*(P(1,3:Nx) - P(1,1:Nx-2))/(2*dx) ...
            );
        W_xx = (2*P(1,1) - 5*P(1,2) + 4*P(1,3) - P(1,4))/dx^2;
        Q(1:1) = (T_0 + Ro*V_0^2) * W_xx ...
            + Ro*V_0 * (-3/2*P(2,1) + 2*P(2,2) - 1/2*P(2,3))/dx;
        Q(2,1) = 0;
        W_xx = (2*P(1,Nx) - 5*P(1,Nx-1) + 4*P(1,Nx-2) - P(1,Nx-3))/dx^2;
        Q(1,Nx) = (T_0 + Ro*V_0^2) * W_xx ...
            + Ro*V_0 * (3/2*P(2,Nx) - 2*P(2,Nx-1) + 1/2*P(2,Nx-2))/dx;
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
                    Q(2,:);
                    zeros(1,Nx)
                ];
            for i = 3:Nx
                M_(:,:,i) = -(A_i_ * L_(:,:,i-1) + B_i_) \ (A_i_ * M_(:,:,i-1) - F(:,i-1));
            end
            Q(:,Nx) = 0;
            Q(3,Nx) = (2 * M_(2,1,Nx)) / (dx^2 - 2*L_(2,3,Nx));
            for i = Nx-1:-1:1
                Q(:,i) = L_(:,:,i+1) * Q(:,i+1) + M_(:,:,i+1);
            end
            Q(1,1) = 0;
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
    elseif k == iterations
        if Ew(k-1) < Ew(k)
            Control = prevControl;
            k = k-1;
        end
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
        FullEnergyTitle(Ew, k);
        SaveAsGif(folder, ['control' image_type], 1, 0);
    end
    k = k+1;
end