close all;
clearvars;
tFormatStr = '%.3f';    %������ ������ ������� �� �������
V0 = 1;                 %���������� �������� �������� ��������
c = sqrt(2);            %C - �������� ��������������� �����
dx = 0.01;              %h - ��� �� ������
dt = 0.0001/2;          %tau - ��� �� �������
Lx = 1;                 %������ ����� ��������
D = 1;                  %? - ����������� ��������
S = Lx;                 %? - ����� ��������(����������)
Ro = 1;                 %��������� ��������� ��������
m = Ro*S;               %����� �� ������� ����� ��������
T0 = 1;                 %���������� ���������� �������� (�� ������� �����)
T = 0.01;               %����� (������� �������)
X = 0:dx:Lx;            %����� �� ����� ��������
Nx = size(X,2);         %����� ����� � ����� X
g_max = 1;              %max ����������� �������
g_max = g_max/m;
g = @(t, val) 0;        %����������� �������

f = NewFigure('������ ������');
saver =  CSaveAsGif();

