classdef CSaveAsGif
    properties (GetAccess = private, SetAccess = private)
        fileName;
        gifDelayTime;
        format;
        fig;
    end
    methods (Access = public)
        function obj = CSaveAsGif(fn, fig, dt)
            obj.fileName = strcat('../gif/',fn);
            obj.gifDelayTime = dt;
            obj.format = 'gif' ;
            obj.fig = fig;
            obj.SaveGif(1);
        end
        function [ ] = SaveGif( obj, create )
            frame = getframe(obj.fig);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256, 'dither');
            if nargin == 1
                create = 0;
            end
            if create == 1
                imwrite(                    ...
                    imind, cm               ...
                    ,obj.fileName, obj.format       ...
                    ,'Loopcount' , inf       ...
                    ,'DelayTime' , obj.gifDelayTime);
            else
                imwrite(                    ...
                    imind, cm               ...
                    ,obj.filename, obj.format       ...
                    ,'WriteMode', 'append'  ...
                    ,'DelayTime', obj.gifDelayTime);
            end
        end
    end
end