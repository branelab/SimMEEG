classdef dipole < handle
% Class Dipole: create a dipole object for dipole fitting
% USAGE:
%   dip = dipole(); 
%       or
%   dip = dipole(num, mri_pos, mri_transform);
%
% Default Properties: 
% num = 1;
% mri_pos         = [nan nan nan]; % positions X, Y, Z of slice positions of MRI slices of size(anatomy.mri.anatomy)
% dip_pos         = [nan nan nan]; % positions X, Y, Z of positions within the leadfield matrix
% dip_ori         (1,2) double {mustBeInRange(dip_ori, 0, 360)} = [0 90];
% dip_on          (1,1) logical {mustBeMember(dip_on, [0 1])} = 1;
% fit_on          (1,1) logical {mustBeMember(fit_on, [0 1])} = 1;
% fit_order       (1,1) double {mustBeInRange(fit_order, 0, 1e10)} = 1;
% sym_type        (1,1) string {mustBeMember(sym_type, ["free", "fixed", "symmetric", "bound"])} = "free";
% sym_idx         (1,1) double {mustBeInRange(sym_idx, 0, 1e10)} = 0;
% sym_dist        (1,1) double {mustBeInRange(sym_dist, -1e10, 1e10)} = 10;
% ori_type        (1,1) string {mustBeMember(ori_type, ["free", "fixed"])} = "fixed";
% fit_ori         (1,1) logical {mustBeMember(fit_ori, [0 1])} = 1;
% vect_idx        (1,1) logical {mustBeMember(vect_idx, [0 1])} = 1;
% search_pos      = [nan nan nan]; % original position X, Y, Z to start search
% search_radius   (1,1) double {mustBeInRange(search_radius, -1e20, 1e20)} = 200;
% dip_clr         (1,3) double {mustBeInRange(dip_clr, 0, 1)} = [0 0 1];
% leadfield       = []; % leadfield matrix [sensors x dip]; 1dip=scalar 3dip=vector
% wts             = []; % pinv(leadfield) [sensor x dip]; 1dip=scalar 3dip=vector
% swf             = [];
% mrk_handles
% wav_handles

    properties  %% Accessible & Changeable
        num = 1;
        mri_pos         = [nan nan nan]; % positions X, Y, Z of slice positions of MRI slices of size(anatomy.mri.anatomy)
        dip_pos         = [nan nan nan]; % positions X, Y, Z of positions within the leadfield matrix
        dip_ori         (1,2) double {mustBeInRange(dip_ori, -360, 360)} = [0 90];
        dip_on          (1,1) logical {mustBeMember(dip_on, [0 1])} = 1;
        fit_on          (1,1) logical {mustBeMember(fit_on, [0 1])} = 1;
        fit_order       (1,1) double {mustBeInRange(fit_order, 0, 1e10)} = 1;
        sym_type        (1,1) string {mustBeMember(sym_type, ["free", "fixed", "symmetric", "bound"])} = "free";
        sym_idx         (1,1) double {mustBeInRange(sym_idx, 0, 1e10)} = 0;
        sym_dist        (1,1) double {mustBeInRange(sym_dist, -1e10, 1e10)} = 10;
        ori_type        (1,1) string {mustBeMember(ori_type, ["free", "fixed"])} = "fixed";
        fit_ori         (1,1) logical {mustBeMember(fit_ori, [0 1])} = 1;
        vect_idx        (1,1) logical {mustBeMember(vect_idx, [0 1])} = 1;
        search_pos      = [nan nan nan]; % original position X, Y, Z to start search
        search_radius   (1,1) double {mustBeInRange(search_radius, -1e20, 1e20)} = 200;
        dip_clr         (1,3) double {mustBeInRange(dip_clr, 0, 1)} = [0 0 1];
        leadfield       = []; % leadfield matrix [sensors x dip]; 1dip=scalar 3dip=vector
        wts             = []; % pinv(leadfield) [sensor x dip]; 1dip=scalar 3dip=vector 
        swf             = [];
        mrk_handles
        ori_handles
        wav_handles
        %% Keep adding properties above to update all current and future "dipole" objects
    end

    methods

        %% Instantiate "dipole" with inputs from user - else it creates default 
        function obj = dipole(num, mri_pos, mri_transform)
            if nargin==3    
                obj.num = num; 
                obj.mri_pos = mri_pos;
                obj = obj.get_dip_pos(mri_transform);
            end
        end

        %% get mri_pos
        function obj = get_mri_pos(obj, mri_transform)
        obj.mri_pos = round( ft_warp_apply(inv(mri_transform), obj.dip_pos) ); 
        
        end
        %% get dip_pos
        function obj = get_dip_pos(obj, mri_transform)
            obj.dip_pos = ft_warp_apply(mri_transform, obj.mri_pos);
        end
        %% calculate source waveform .swf
        function obj = calc_swf(obj,data)
            for v=1:size(obj,2) % loop through dipoles
                if ~isempty(obj(v).wts)
                    obj(v).swf = [];
                    for t=1:size(data,3) % trials --> data=[samp x chans x trials];
                        obj(v).swf(:,:,t) = obj(v).wts'*squeeze(data(:,:,t))';
                    end
                end

            end
        end

    end
end