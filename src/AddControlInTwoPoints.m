function [ g ] = AddControlInTwoPoints(val, i, abs_val)
    g = 0;
    if i==1
        if val > 0
            g = abs_val;
        end
    elseif i==2
        if val < 0
            g = -abs_val;
        end
    end
end