classdef CSaveAsGif
    properties (GetAccess = public, SetAccess = private)
        fileName;
        gifDelayTime;
        format;
    end
    methods (Access = public)
        function obj = CSaveAsGif(fn, dt)
            obj.fileName = strcat('../gif/',fn);
            obj.gifDelayTime = dt;
            obj.format = 'gif' ;
            obj.SaveGif(1);
        end
        function [ ] = SaveGif( obj, create )
            frame = getframe(1);
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