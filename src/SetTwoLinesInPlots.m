function [ ] = SetTwoLinesInPlots( p_1, p_2, data, t, name_1, name_2)
    p_1.YData = data(1,:);
    p_2.YData = data(2,:);
    Legend(t, name_1, name_2);
    drawnow;
end

