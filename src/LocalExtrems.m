function [Min, Max] = LocalExtrems( Values, Min, Max )
    LocalMin = min(Values(1,:));
    LocalMax = max(Values(1,:));
    if LocalMin < Min
        Min = LocalMin;
    end
    if LocalMax > Max
        Max = LocalMax;
    end
    global FormatStr;
    title({
        ['Min (W) = ' num2str(Min,FormatStr)]
        ['LocalMin (W) = ' num2str(LocalMin,FormatStr)]
        ['Max (W) = ' num2str(Max,FormatStr)]
        ['LocalMax (W) = ' num2str(LocalMax,FormatStr)]})
end

