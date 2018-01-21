function [ ] = Legend( t, f, s )
    global FormatStr;
    switch nargin
        case 1
            legend(['t = ' num2str(t, FormatStr)]);
        case 2
            legend(['t = ' num2str(t, FormatStr)], f);
        case 3
            legend(['t = ' num2str(t, FormatStr)], f, s);
    end
end

