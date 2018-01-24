function [ ] = SetTwoLinesInPlots( p_1, p_2, data, t, name_1, name_2 )
    switch  nargin 
        case 5
            p_1.YData = data(2,:);
            Legend(t, name_1);
        case 6
            p_1.YData = data(1,:);
            p_2.YData = data(2,:);
            Legend(t, name_1, name_2);
        otherwise    %для отрисовки управления в точках
            scatter(data, p_1, 'o', 'MarkerEdgeColor','r', 'LineWidth',1);
            scatter(data, p_2, '+', 'MarkerEdgeColor','b', 'LineWidth',1);
            legend('g_1', 'g_2');
    end
    drawnow;
end

