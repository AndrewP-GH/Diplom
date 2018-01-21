classdef CSaveAsGif
    properties(Access = private)
        fileName = 'task';
        gifDelayTime = 0;
        format = 'gif';
        path = '../gif/';
    end
    properties
        F_ind = 1;
    end
    methods
        function obj = CSaveAsGif(fn, f_ind, dt)
            if nargin == 1
                obj.fileName = obj.GetFullFileName(fn);
            else
                obj.fileName = obj.GetFullFileName(obj.fileName);
            end
            if nargin == 2
                obj.F_ind = f_ind;
            end
            if nargin == 3
                obj.gifDelayTime = dt;
            end
            obj.SaveGif(1);
        end
        function [ ] = SaveGif( obj, create )
            [imind, cm] = CSaveAsGif.GetImg(obj.F_ind);
            if nargin == 1
                create = 0;
            end
            if create == 1
                imwrite(                        ...
                    imind, cm                   ...
                    ,obj.fileName, obj.format   ...
                    ,'WriteMode', 'overwrite'   ...
                    ,'Loopcount' , Inf          ...
                    ,'DelayTime' , obj.gifDelayTime);
            else
                imwrite(                        ...
                    imind, cm                   ...
                    ,obj.fileName, obj.format   ...
                    ,'WriteMode', 'append'      ...
                    ,'DelayTime', obj.gifDelayTime);
            end
        end
        function [ ] = SaveTif( obj, fName , f_ind)
            imFormat = 'tif';
            if nargin == 2
                f_ind = obj.F_ind;
            end
            [imind, cm] = CSaveAsGif.GetImg(f_ind);
            fName = obj.GetFullFileName(fName, imFormat);
            imwrite(imind, cm                   ...
                    ,fName, imFormat            ...
                    ,'WriteMode', 'overwrite');
        end
        function fName = GetFullFileName( obj, name, format )
            if nargin == 2
                format = obj.format;
            end
            fName = strcat(obj.path, name, '.', format);
        end
    end
    methods(Static)
        function [imind, cm] = GetImg(f_ind)
            frame = getframe(f_ind);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256, 'dither');
        end
    end
end