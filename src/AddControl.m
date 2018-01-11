function [ F ] = AddControl(F, val, g)
    global dx;
    global dt;
    g = 0;
    if val > 0
        g = 0.5;
    elseif val < 0
        g = -0.5;
    end
    if g~=0
        F = F + [0; dx^2*dt*g];
    end
end

