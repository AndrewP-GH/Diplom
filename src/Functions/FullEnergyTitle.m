function [] = FullEnergyTitle( Ew, k )
    t = title({
            strcat('Ew_{νΰχ.} =', [' ' num2str(Ew(1))])
            strcat('Ew_{ξος.} =', [' ' num2str(Ew(k))])
        });
    set(t, 'horizontalAlignment', 'left')
end

