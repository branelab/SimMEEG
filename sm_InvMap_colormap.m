function sm_InvMap_colormap(varargin)
% set colormap for 3D InvMap axis
global h
colormap(h.axes_3D_images,h.menu_InvMap_colormap.String{h.menu_InvMap_colormap.Value}); 