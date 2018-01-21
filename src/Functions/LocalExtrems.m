function [Min, Max] = LocalExtrems( Values, Min, Max, Raw, Name )
    LocalMin = min(Values(Raw,:));
    LocalMax = max(Values(Raw,:));
    if LocalMin < Min
        Min = LocalMin;
    end
    if LocalMax > Max
        Max = LocalMax;
    end
    global FormatStr;
    title({
        ['Min (' Name ') = ' num2str(Min,FormatStr)]
        ['LocalMin (' Name ') = ' num2str(LocalMin,FormatStr)]
        ['Max (' Name ') = ' num2str(Max,FormatStr)]
        ['LocalMax (' Name ') = ' num2str(LocalMax,FormatStr)]})
end

