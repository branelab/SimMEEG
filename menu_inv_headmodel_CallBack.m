
function menu_inv_headmodel_CallBack(varargin)
global h

if h.menu_inv_headmodel.Value == 1     % Whole brain
%     h.anatomy.leadfield_inv_model = h.anatomy.vol_leadfield;
        h.menu_inv_ori_normal.Enable = 'inactive';  h.menu_inv_ori_normal.ForegroundColor=[1 1 1]*.5;   h.menu_inv_ori_normal.Value = 1;    % Random orientations locked
elseif h.menu_inv_headmodel.Value == 2     % Cortical Surface
%     h.anatomy.leadfield_inv_model = h.anatomy.mesh_leadfield;
%     h.menu_inv_ori_normal.Enable = 'on';  h.menu_inv_ori_normal.ForegroundColor=[1 1 1]*0; h.menu_inv_ori_normal.Value = 2;    % Random orientations locked
    h.menu_inv_ori_normal.Enable = 'inactive';  h.menu_inv_ori_normal.ForegroundColor=[1 1 1]*.5; h.menu_inv_ori_normal.Value = 1;    % Random orientations locked
end