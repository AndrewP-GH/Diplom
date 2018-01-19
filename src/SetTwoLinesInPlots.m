function [ ] = SetTwoLinesInPlots( p_1, p_2, data, t, name_1, name_2 )
    if nargin == 6
        p_1.YData = data(1,:);
        p_2.YData = data(2,:);
        Legend(t, name_1, name_2);
    else
        scatter(data, p_1, 'o', 'MarkerEdgeColor','r', 'LineWidth',1);
        scatter(data, p_2, '+', 'MarkerEdgeColor','b', 'LineWidth',1);
        legend('g_1', 'g_2');
    end
    drawnow;
end

