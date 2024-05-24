function sm_ARM_change_source_params(varargin)
global h

switch varargin{end}
    case 'Open'
        if  h.btn_ARM_change_source_params.Value == 1
            h.panel_ARM_source_params.Visible = 'on';
        else
            h.panel_ARM_source_params.Visible = 'off';
        end
    case 'Update'
        % convert Orientations to Az and El
        [Az,El] = cart2sph(h.cfg.ARM_params.vx_ori(:,1),h.cfg.ARM_params.vx_ori(:,2),h.cfg.ARM_params.vx_ori(:,3));
        OriAz = rad2deg(Az); OriEl = rad2deg(El);
        OriAz(OriAz<0) = 360+OriAz(OriAz<0); % rescaling to be between 0 - 360 degrees
        OriEl(OriEl<0) = 360+OriEl(OriEl<0); % rescaling to be between 0 - 360 degrees
        
        data = [h.cfg.ARM_params.vx_idx' round(h.cfg.ARM_params.vx_locs) round(OriAz) round(OriEl) h.cfg.ARM_params.vx_amp' h.cfg.ARM_params.sig_amp_perc' h.cfg.ARM_params.prepost_amp_perc' h.cfg.ARM_params.sig_latency h.cfg.ARM_params.sig_risetime];
        h.table_ARM_source_params.Data = data;
    case 'Apply'
        vx_locs = h.table_ARM_source_params.Data(:,2:4);
        h.cfg.ARM_params.vx_idx = find_nearest_voxel(vx_locs,h.anatomy.leadfield.voxel_pos);
        h.table_ARM_source_params.Data(:,1) = h.cfg.ARM_params.vx_idx;
        h.table_ARM_source_params.Data(:,2:4) = round(h.anatomy.leadfield.voxel_pos(h.cfg.ARM_params.vx_idx,1:3));
        data = h.table_ARM_source_params.Data;
        h.cfg.ARM_params.vx_idx = data(:,1)';
        h.cfg.ARM_params.vx_locs = data(:,2:4);
        [x,y,z] = sph2cart(deg2rad(data(:,5)),deg2rad(data(:,6)),ones(size(data(:,5))));
        h.cfg.ARM_params.vx_ori = [x y z];
        h.cfg.ARM_params.vx_amp = data(:,7)';
        h.cfg.ARM_params.sig_amp_perc = data(:,8)';
        h.cfg.ARM_params.prepost_amp_perc = data(:,9)';
        h.cfg.ARM_params.sig_latency = data(:,10:11);
        h.cfg.ARM_params.sig_risetime = data(:,12);
        h.fcn_handle.plot_3D_mri('');
    case 'Load Excel'
        answ = questdlg(sprintf('Loading will overwrite the ARM Source Parameters.\n\nDo you want to overwrite?\n'),'Overwrite Parameters?','Yes','No','No');
        switch answ
            case 'Yes'
                if ~isfield(h.cfg,'ARM_params.saved_params_file')
                    h.cfg.ARM_params.saved_params_file = fullfile(h.data_dir, sprintf('%s_ARM_params.xlsx',h.cfg.study.study_name));
                end
                [fname,fpath] = uigetfile({'*.xlsx'},'Save ARM Source Parameters',h.cfg.ARM_params.saved_params_file); % user can change name and path to save
                h.cfg.ARM_params.saved_params_file = fullfile(fpath,fname);
                
                h.waitfor_panel.Visible = 'on';
                h.waitfor_txt.String = sprintf('Saving File: %s\n',h.cfg.ARM_params.saved_params_file);
                
                [data] = xlsread(h.cfg.ARM_params.saved_params_file);
                if size(data,2)==12
                    h.table_ARM_source_params.Data = data;
                else
                    txt = sprintf('   %s\n',h.table_ARM_source_params.ColumnName{:});
                    warndlg(sprintf('ARM Params Not Chanaged!\nFile does not have correct format\nRows = Sources\nColumns = \n%s\n',txt),'');
                end
                
                h.waitfor_txt.String = 'Default Message';
                h.waitfor_panel.Visible = 'off';
                
            case 'No'
        end
    case 'Save Excel'
        if ~isfield(h.cfg,'ARM_params.saved_params_file')
            h.cfg.ARM_params.saved_params_file = fullfile(h.data_dir, sprintf('%s_ARM_params.xlsx',h.cfg.study.study_name));
        end
        [fname,fpath] = uiputfile({'*.xlsx'},'Save ARM Source Parameters',h.cfg.ARM_params.saved_params_file); % user can change name and path to save
        h.cfg.ARM_params.saved_params_file = fullfile(fpath,fname);
        
        h.waitfor_panel.Visible = 'on';
        h.waitfor_txt.String = sprintf('Saving File: %s\n',h.cfg.ARM_params.saved_params_file);
        
        clmn_names = h.table_ARM_source_params.ColumnName;
        data = h.table_ARM_source_params.Data;
        c(1,:) = clmn_names;
        c(2:1+size(data,1),:) = num2cell(data);
        xlswrite(h.cfg.ARM_params.saved_params_file,c);
        h.saved_ARM_source_params_txt.String = sprintf('Saved ARM Parameter File: %s',fname);
        
        h.waitfor_txt.String = 'Default Message';
        h.waitfor_panel.Visible = 'off';
       
end


