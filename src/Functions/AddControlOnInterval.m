function [ Control ] = AddControlOnInterval( intervals, X, QH, g_max)
    Control = zeros(1, length(X));
    n = length(intervals);
    for i = 1:n
        int = cell2mat(intervals(i));
        integral = trapz(X(int), QH(int));
        m = mod(i,2);
        if m == 1 && integral > 0
            Control(int) = g_max;
        elseif m == 0 && integral < 0
            Control(int) = -g_max;
        end 
    end
end

