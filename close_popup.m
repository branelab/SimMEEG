function close_popup(~,~,h_fig)
global h
switch h_fig.Name
    case 'Brainstorm Anatomy'
        m = findobj(h_fig.Children,'Style','popupmenu');
        h.popup_data = m.String{m.Value};
    case ''
end
close(h_fig)

