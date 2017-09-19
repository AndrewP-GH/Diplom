function [ output_args ] = mainSystem( dx, dt, D, m, W, V, U, T, Lx )
    P = [ W; V; U];
    P(:,1) = 0;                         %граничное левое
    P(:,Nx) = 0;                        %граничное правое

    Alpha = zeros(3,3,Nx);
    Beta = zeros(3,1,Nx);

    nu = dt/(dx^2);
    A = [ 
            0 0 0; 
            nu*(V0^2-c^2) -V0*dt/dx D/m*nu; 
            1 0 0 
        ];
    B = [ 
            1 -dt 0; 
            -2*nu*(V0^2-c^2) 1 -2*D/m*nu; 
            -2 0 -dx^2 
        ];
    C = [ 
            0 0 0; 
            nu*(V0^2-c^2) V0*dt/dx D/m*nu; 
            1 0 0 
        ];

    hold all;
    axis([0 Lx 0 0.2]);
    grid on;
    plot(0,0, 'w');
    p1 = plot(X,P(1,:));
    set(p1,'LineWidth', 4); 
    legend(['t=' num2str(0,tFormatStr)],'W');
    set (p1, 'Xdata', X, 'Ydata', P(1,:));
    drawnow;  
   
    for t=0:dt:T-dt
        F = [ P(1,:); P(2,:) + g(t,val_g)*dt; zeros(1,Nx) ];
        Alpha(:,:,2) = [ 0 0 0; 0 0 0; 0 2/(dx^2) 0];
        for i=3:Nx
            Alpha(:,:,i) = -(A*Alpha(:,:,i-1)+B)\C;
            Beta(:,:,i) = (A*Alpha(:,:,i-1)+B)\(F(:,i-1)-A*Beta(:,:,i-1));
        end
        for i=Nx-1:-1:2
            P(:,i) = Alpha(:,:,i+1)*P(:,i+1)+Beta(:,:,i+1);
        end
        set (p1, 'Xdata', X, 'Ydata', P(1,:));
        legend(['t=' num2str(t+dt,tFormatStr)],'W');

        drawnow;  
        SaveAsGif(filename, gifDelayTime, 1);
    end
end

