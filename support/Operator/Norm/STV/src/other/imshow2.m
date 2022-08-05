function imshow2(u)
% IMSHOW2 (display image into a Matlab Figure) with tight borders AND without
% any interpolation (one pixel of the screen = one pixel of the image),
% which is unfortunately not the case of imshow.

    function callback_keyboard_pressed(fig_hdl,evt) % keyboard event handler function
        if strcmp(evt.Key,'q') % stop process & quit
            close(fig_hdl);
        end
    end

    gr = groot();
    u = double(u);
    [ny,nx] = size(u);
    image(u,'CDataMapping','scaled','XData',0:nx-1,'YData',0:ny-1);
    axis 'image';
    fg = gcf();
    set(fg,'Colormap',gray(256),'KeyPressFcn',@callback_keyboard_pressed); 
    set(gca(),'XLim',[0,nx-1],'YLim',[0,ny-1],'Visible','off','Units','pixels','ActivePositionProperty','position','Position',[1,1,nx,ny]);
    fg.Position(1:2) = [min(100,gr.ScreenSize(3)-nx),min(gr.ScreenSize(4)-ny-100,gr.ScreenSize(4)-ny)];
    fg.Visible = 'on';
    pause(1e-6);
    fg.Position(3:4) = [nx,ny];

end
