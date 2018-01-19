function [ Ew ] = PaperFullEnergy( P, T_0, Ro, V_0)
    global dx;
    Nx = length(P(1,:));
    Wx = zeros(1,Nx);
    Wx(2:Nx-1) = (P(1,3:Nx)-P(1,1:Nx-2))/(2*dx);
    Wx(1) = (-3/2*P(1,1)+2*P(1,2)-1/2*P(1,3))/dx;
    Wx(Nx) = (3/2*P(1,1)-2*P(1,2)+1/2*P(1,3))/dx;
    Wx = 1/2 * (T_0*Wx.^2 + Ro*(P(2,:) + V_0*Wx).^2);
    Ew = (sum(Wx)-1/2*(Wx(1)+Wx(Nx)))*dx;
end