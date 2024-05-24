function sm_mousemove(hobject,eventdata)
% Monitors mouse moves over different panels, axes, ... and executs commands

global h

cpt = get(h.main_fig,'currentpoint');
display(cpt);


%% Axes: Source waves --> Turn off rotate3D and display point locations
if cpt(1)>h.axes_3D_images.Position(1) && cpt(1)<sum(h.axes_3D_images.Position([1 3])) ...
        && cpt(2)>h.axes_3D_images.Position(2) && cpt(2)<sum(h.axes_3D_images.Position([2 4]))
    rotate3d(h.axes_3D_images,'on')
end


%% Axes: Source waves --> Turn off rotate3D and display point locations
if cpt(1)>h.axes_source_waves.Position(1) && cpt(1)<sum(h.axes_source_waves.Position([1 3])) ...
        && cpt(2)>h.axes_source_waves.Position(2) && cpt(2)<sum(h.axes_source_waves.Position([2 4]))
    if ~isvalid(h.text_cursor_specs)
        h.text_cursor_specs = text(h.axes_source_waves, h.axes_source_waves.XLim(1)+(.02*range(h.axes_source_waves.XLim)) , ...
            h.axes_source_waves.YLim(2)-(.1*range(h.axes_source_waves.YLim)), 'X: Y:','Color','r');
    end
    rotate3d(h.axes_source_waves,'off')
    h.text_cursor_specs.String = sprintf('%.3f, %.3f', cpt(1,1:2));
end



end


