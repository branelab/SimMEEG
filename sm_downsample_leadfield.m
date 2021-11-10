function sm_downsample_leadfield(varargin)
global h

%% spatially down sampling lead field grid by a factor of "ds_factor" times 

ds_factor = str2num(h.edit_inv_plv_lf_grid_idx.String);
if ~(mod(ds_factor,1)==0) || ds_factor<0 % then it's an whole number
    ds_factor = round(abs(ds_factor)); 
    fprintf('Down Sampling factor must be a whole number. Converted it to %.f and downsampled lead field grid.\n',ds_factor);
    h.edit_inv_plv_lf_grid_idx.String = ds_factor;
end



if ds_factor<1
    hm = warndlg('Down Sample factor must be > 1','Please Change');
elseif ds_factor ==1 % reset to full leadfield grid
    h.inv_soln(h.current_inv_soln).plv_leadfield_grid_idx = 1:size(h.inv_soln(h.current_inv_soln).leadfield.voxel_pos,1);
elseif ds_factor>1
    switch h.inv_soln(h.current_inv_soln).headmodel_type
        case 'Whole Brain'
            
            %% creating grid "vx_grid" as per Alex's description "nVoxels = nX * nY * nZ, and the ordering is such that the last index (i.e. Z direction) changes most fast, and the first (i.e. X) - most slow"
            vx_pos = h.inv_soln(h.current_inv_soln).leadfield.voxel_pos;
            vx_res = diff(vx_pos(:,1)); vx_res = min(abs(vx_res(vx_res~=0)));   % voxel resolution
            vx_res = vx_res * ds_factor; % downsampling
            % making 3D grid
            xg = min(vx_pos(:,1))-vx_res:vx_res:max(vx_pos(:,1))+vx_res;
            yg = min(vx_pos(:,2))-vx_res:vx_res:max(vx_pos(:,2))+vx_res;
            zg = min(vx_pos(:,3))-vx_res:vx_res:max(vx_pos(:,3))+vx_res;
%             vx_grid = meshgrid(xg,yg,zg);
            vx_grid=[]; v=0;    
            for d1=1:length(xg)
                for d2=1:length(yg)
                    for d3=1:length(zg)
                        v=v+1;
                        vx_grid(v,:) = [xg(d1) yg(d2) zg(d3)];
                    end
                end
            end
            
            %% find equal to leadfield to new down sampled vx_grid to only include the vx_grid near to the leadfields when searching for nearest leadfield idx
            vx_grid = round(vx_grid); vx_pos=round(vx_pos); % need becuase of decimal error when finding exact locations
            v_idx = []; vv=[];
            for v=1:size(vx_grid,1)
                found_idx = find((vx_grid(v,1)==vx_pos(:,1))==1 & (vx_grid(v,2)==vx_pos(:,2))==1 & (vx_grid(v,3)==vx_pos(:,3))==1);
              if ~isempty(found_idx )
                  v_idx(v)=found_idx; 
              end
            end
            v_idx = v_idx(v_idx~=0); % leadfidel indices
            
            %% finding indices of vx_grid within lead field grid
            h.inv_soln(h.current_inv_soln).plv_leadfield_grid_idx = v_idx;
            
        case 'Cortical Surface'
            
            fac = h.inv_soln(h.current_inv_soln).headmodel_mesh.tri;
            verts2 = h.inv_soln(h.current_inv_soln).headmodel_mesh.pos;
            hf = figure(4); clf; scalp2 = patch('Faces',fac,'Vertices',verts2, 'FaceColor',[1,.75,.65],'EdgeColor','none'); view(135,30); axis tight; daspect([1,1,1]); lightangle(145,30); lighting phong;
            close(hf); 
            fv.vertices=verts2; fv.faces=fac; ds_mesh=(size(verts2,1)/ds_factor); nfv = reducepatch(fv,ds_mesh); % downsample faces
            mesh_volumes(1).pos=nfv.vertices; mesh_volumes(1).tri=nfv.faces;
            v_idx = find_nearest_voxel(nfv.vertices,h.inv_soln(h.current_inv_soln).leadfield.voxel_pos);
            h.inv_soln(h.current_inv_soln).plv_leadfield_grid_idx = v_idx;
            
            
    end
    
end

sm_menu_inv_analyses_CallBack; % updating comparison locs in case "Lead Field Grid" is selected


if exist('hm','var'); delete(hm); end

