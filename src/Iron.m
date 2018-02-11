clearvars; clc; close all;
addpath('./Functions');
Figure = NewFigure('Железо');

%% Общие параметры задачи
global FormatStr        %формат вывода времени на графике
FormatStr = '%.3f';
V_0 = 0.5;              %продольная скорость движения пластины
C = sqrt(2);            %C - скорость распространения волны
V_C = V_0^2 - C^2;      
global dx;
global dt;
dx = 0.01;              %h - шаг по длинне
dt = 0.01;              %tau - шаг по времени
L_x = 1;                %длинна между роликами
D = 1;                  %? - неизвестный параметр
S = L_x;                %? - длина пластины(абсолютная?)
Ro = 1;                 %плотность материала пластины
m = Ro*S;               %масса на единицу длины пластины
T_0 = 1;                %продольное напряжение пластины (на единицу длины)

T = 1;                  %время (верхняя граница)
Tarray = 0:dt:T;
Nt = T/dt+1;            %число слоев повремени

X = 0:dx:L_x;           %сетка по длине пластины
Nx = size(X, 2);        %число узлов в сетке X
g_max = -3000;              %max управляющей функции
i_1 = round(Nx/3);
i_2 = round(Nx/3*2);
iterations = 2;         %число итераций
Ymax = 2;

global CalcExtrems      %выводить минимум и максимум
CalcExtrems = true;
OnlyQ2 = true;

folder = CreateImageFolder([datestr(now, 'dd-mmm-yyyy HH_MM_SS') '_железо']);
image_type = '.tiff';
gif_type = '.gif';
image_name = 'tmp';
save_gif = true;
gif_delay = 1/24;

%% Параметры прямой задачи
W = 0.2*sin(pi*X);	    %начальное распределение точек листа
V = zeros(1,Nx);        %начальное распределение скорости листа
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
L = zeros(3,3,Nx);      %прогоночные коэффициенты прямой задачи
M = zeros(3,1,Nx);      %прогоночные коэффициенты прямой задачи

L(:,:,3) = - B_i \ C_i;
for i = 4:Nx
    L(:,:,i) = -(A_i * L(:,:,i-1) + B_i) \ C_i;
end

%% Параметры обратной задачи
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
L_ = zeros(3,3,Nx);     %прогоночные коэффициенты обратной задачи
M_ = zeros(3,1,Nx);     %прогоночные коэффициенты обратной задачи
Q = zeros(3,Nx);        %вектор значений (q_1, q_2, q_3)

L_(:,:,2) = [   0 0      0; 
                0 0      0; 
                0 2/dx^2 0 
            ];
for i = 3:Nx
    L_(:,:,i) = -(A_i_ * L_(:,:,i-1) + B_i_) \ C_i_;
end

QH = zeros(Nt,Nx);      %тут будут храниться значения q2 в каждый момент времени
Control = zeros(2,Nt);

%% Цикл вычисления
for k=1:iterations
    %% Вычисление прямой задачи
    [p_1, p_2] = FigurePrepare(Figure, X);
    axis([0 L_x -Ymax Ymax]);
    p_1.LineWidth = 3;
    P = [W; V; U];         %вектор значений (W, V)
    P(:,1) = 0;             %граничное левое
    P(:,Nx) = 0;            %граничное правое
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
        Control(1,t+1) = AddControlInTwoPoints(QH(t,i_1), 1, g_max);
        Control(2,t+1) = AddControlInTwoPoints(QH(t,i_2), 2, g_max);
        F(2,i_1) = F(2,i_1) + dt*Control(1,t+1);
        F(2,i_2) = F(2,i_2) + dt*Control(2,t+1);
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
    disp([ 'Full energy Ew = ' num2str(PaperFullEnergy(P, T_0, Ro, V_0)) ]);
    SaveAsGif(folder, [image_name '_end' image_type], 1, 0);
    %% Вычисление обратной задачи
    if k ~= iterations
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
    else
        FigurePrepare(Figure, Tarray, false);
        SetTwoLinesInPlots(Control(1,:), Control(2,:), Tarray);
        SaveAsGif(folder, ['control' image_type], 1, 0);
    end
end