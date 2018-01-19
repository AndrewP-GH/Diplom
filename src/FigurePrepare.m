function [p_1, p_2] = FigurePrepare(Figure, X, draw_line)
    clf(Figure);
    hold all;
    grid on;
    if nargin ~= 3 || draw_line == true
        plot(0,0, 'w');         %��� ����������� ������� � �������
        p_1 = plot(0,0);        %���������� ����� 1
        p_2 = plot(0,0);        %���������� ����� 2
        p_1.XData = X;
        p_2.XData = X;
    end
    global CalcExtrems
    if CalcExtrems == true
        title({
        ' '
        ' '
        ' '
        ' '});
    end
end

