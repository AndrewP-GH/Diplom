function [] = FullEnergyTitle( Ew, k )
    t = title({
            strcat('Ew_{���.} =', [' ' num2str(Ew(1))])
            strcat('Ew_{���.} =', [' ' num2str(Ew(k))])
        });
    set(t, 'horizontalAlignment', 'left')
end

