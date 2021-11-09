function set_3D_transparency(varargin)
global h
alpha_gain = [.25 1 .5];
switch h.inv_soln(h.current_inv_soln).headmodel_type
    case 'Whole Brain' % whole brain
        xh=handle(h.func3D);
        for a=1:length(xh); xh(a).FaceAlpha = h.slider_3D_transparency_func.Value; end
        for a=1:length(h.anat3D); h.anat3D(a).FaceAlpha = h.slider_3D_transparency_anat.Value*alpha_gain(a); end
    case 'Cortical Surface'    % Cortical Surface
        h.func3D(1).FaceAlpha = h.slider_3D_transparency_anat.Value*alpha_gain(1);
        h.func3D(2).FaceAlpha = h.slider_3D_transparency_anat.Value*alpha_gain(2);
        h.func3D(3).FaceAlpha = h.slider_3D_transparency_func.Value; 
end
