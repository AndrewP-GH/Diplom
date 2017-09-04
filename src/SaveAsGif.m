function [ ] = SaveAsGif( filename, gifDelayTime, add )
    filename = strcat('../gif/',filename);
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256, 'dither');
    format = 'gif' ;
    if add == 0                     %create new gif
        imwrite(                    ...
            imind, cm               ...
            ,filename, format       ...
            ,'Loopcount' ,inf       ...
            ,'DelayTime' ,gifDelayTime);
    else                            %append gif
        imwrite(                    ...
            imind, cm               ...
            ,filename, format       ...
            ,'WriteMode', 'append'  ...
            ,'DelayTime', gifDelayTime);
    end
end
