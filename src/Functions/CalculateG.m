function [ g_min, g, g_max ] = CalculateG( k, h_g, g_min, g, g_max, Ew )
    switch k    %����� ����������� ��������
        case 1  %��� ����������
            g_min = g;
            g = g_max;
        case 2  %
            g_max = g_max + h_g;
            g_min = g;
            g = g_max;
        case 3
            
        otherwise
            
    end
    
end

