function [p_1, p_2] = FigurePrepare(Figure, X)
    clf(Figure);
    hold all;
    grid on;
    plot(0,0, 'w');         %для отображения времени в легенде
    p_1 = plot(0,0);        %дескриптор линии 1
    p_2 = plot(0,0);        %дескриптор линии 2
    p_1.XData = X;
    p_2.XData = X;
    global CalcExtrems
    if CalcExtrems == true
        title({
        ' '
        ' '
        ' '
        ' '});
    end
end

