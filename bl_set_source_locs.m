function h2=bl_set_source_locs(mri,source_locs,vol)
%   This program will plot mri and allow for selecting 3 sources to be entered into SimSignals.m software.
%   Source locations can either be set by moving the cursors to a specific location and then clicking on "Set Source 1"
%   or by clicking and dragging any of the 3 sources. Remember to click on "Save to Workspace" once the 3 sources are in 
%   the locations you want to simulate. 
%
% INPUT
%   mri = Field Trip's structure of loading mri. Note, that the transformation matrix mri.transform should be the same one used for creating the headmodel and leadfields.
%           must contain:
%               mri.dim = dimensions of mri.anatomy
%               mri.hdr = header information read in using Field Trip
%               mri.anatomy = in voxel space coordinates
%               mri.transform = 4x4 tranformation matrix to convert voxel to mm coordinates.
%   source_locs = location in mri.anatomy space coordinates
%   vol = (optional) volume data for source model to plot source_locs_mm in 3D (see bl_plot_source_locs.m for more options)
%
% OUTPUT
%   Click on "Save to Workspace" to load the source locations into matlab's base workspace. 
%       source_locs = source locations in voxel coordinates
%       source_locs_mm = source locations in mm after transformation using ft_warp_apply.m (see Field Trip for more information)
%
%% initialize figure
global h2
h2=figure(1); h2.Color='w'; h2.Position=[300 100 850 700];


if nargin<2 || isempty(source_locs)
    h2.UserData.source_locs=[35 109 78; 151 109 78; 91 48 78]; %ones(3,3); %[];
else
    h2.UserData.source_locs=source_locs;
end
h2.UserData.x_sagittal_slice=h2.UserData.source_locs(1,1); % set starting cursor to first source loc
h2.UserData.y_coronal_slice=h2.UserData.source_locs(1,2);
h2.UserData.z_axial_slice=h2.UserData.source_locs(1,3);

% h2.UserData.fiducial=mri.hdr.fiducial;  % Not implementing plotting of fiducials yet
h2.UserData.transform=mri.transform;
h2.UserData.source_clr=[0 0 1; 1 0 0; 0 0.8 0];
h2.UserData.source_size=60;
h2.UserData.anatomy=mri.anatomy;
h2.UserData.dims=size(h2.UserData.anatomy);
h2.UserData.drag_axial_source=0; h2.UserData.drag_sagittal_source=0; h2.UserData.drag_coronal_source=0;

%% Panel 1: MRI information - Upper Left Panel
hp1 = uipanel(h2,'Title','MRI Information','FontSize',12,...
    'BackgroundColor','w','Foregroundcolor','k',...
    'Position',[.01 .51 .47 .47],'Units','normalize');
mri_txt=uicontrol(hp1,'Style','text',...
    'BackgroundColor','w','Foregroundcolor','k',...
    'Position',[.01 .51 .45 .45],'Units','normalize',...
    'String',''); mri_txt.Position=[.01 .81 .9 .15];
h2.UserData.mri_str=sprintf('MRI File: %s\n',mri.hdr.mrifile);
mri_txt.String=h2.UserData.mri_str; mri_txt.HorizontalAlignment='left';
%% text output
loc_txt=uicontrol(hp1,'Style','text',...
    'BackgroundColor','w','Foregroundcolor','k',...
    'Position',[.01 .61 .55 .55],'Units','normalize',...
    'String','');loc_txt.Position=[.01 .31 .9 .55];
xyz=[h2.UserData.x_sagittal_slice,h2.UserData.y_coronal_slice,h2.UserData.z_axial_slice];
xyz_w=ft_warp_apply(h2.UserData.transform,[1 1 1]); h2.UserData.current_loc=xyz_w;
h2.UserData.source_locs_mm=ft_warp_apply(h2.UserData.transform,h2.UserData.source_locs); h2.UserData.current_loc=xyz_w;
h2.UserData.loc_str=sprintf('pixel =   %.f   %.f   %.f\n mm =   %.f   %.f   %.f\n\nSource 1 =    %.f   %.f   %.f\nSource 2 =    %.f   %.f   %.f\nSource 3 =    %.f   %.f   %.f\n\nSource 1 =    %.f   %.f   %.f\nSource 2 =    %.f   %.f   %.f\nSource 3 =    %.f   %.f   %.f',...
    xyz,xyz_w,h2.UserData.source_locs',h2.UserData.source_locs_mm');
loc_txt.String=h2.UserData.loc_str; loc_txt.HorizontalAlignment='left';
%% btn for source 1 selection locations
btn1 = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(1,:),'ForegroundColor',[1 1 1]*.6,'Style','pushbutton','String','Set Source 1','Position',[.15 .05 .15 .25],'units','normalize','Callback',@select_location); btn1.Position=[.02 .23 .41 .1];
btn2 = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(2,:),'ForegroundColor',[1 1 1]*0,'Style','pushbutton','String','Set Source 2','Position',[.45 .05 .15 .25],'units','normalize','Callback',@select_location); btn2.Position=[.02 .12 .41 .1];
btn3 = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(3,:),'ForegroundColor',[1 1 1]*0,'Style','pushbutton','String','Set Source 3','Position',[.75 .05 .15 .25],'units','normalize','Callback',@select_location); btn3.Position=[.02 .01 .41 .1];
%% btn to go to source 1 selection locations
btn1_loc = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(1,:),'ForegroundColor',[1 1 1]*.6,'Style','pushbutton','String','Go to Source 1','Position',[.15 .05 .15 .25],'units','normalize','Callback',@goto_location); btn1_loc.Position=[.52 .23 .41 .1];
btn2_loc = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(2,:),'ForegroundColor',[1 1 1]*0,'Style','pushbutton','String','Go to Source 2','Position',[.45 .05 .15 .25],'units','normalize','Callback',@goto_location); btn2_loc.Position=[.52 .12 .41 .1];
btn3_loc = uicontrol(hp1,'BackgroundColor',h2.UserData.source_clr(3,:),'ForegroundColor',[1 1 1]*0,'Style','pushbutton','String','Go to Source 3','Position',[.75 .05 .15 .25],'units','normalize','Callback',@goto_location); btn3_loc.Position=[.52 .01 .41 .1];
%% btn save source locs to workspace
btn_save = uicontrol(hp1,'BackgroundColor',[1 .6 0],'Style','pushbutton','String','Save Locations','units','normalize','Callback',@save_location); btn_save.Position=[.652 .88 .3 .1];
%% btn default source locs to workspace
btn_defaults = uicontrol(hp1,'BackgroundColor',[1 1 1]*.8,'Style','pushbutton','String','Default Locs','units','normalize','Callback',@default_location); btn_defaults.Position=[.652 .68 .3 .1];
%% btn for plotting source locs in 3D volume
btn_3D = uicontrol(hp1,'BackgroundColor',[1 .6 .6],'Style','pushbutton','String','Plot 3D','units','normalize','Callback',@plot_source_locs_3D); btn_3D.Position=[.652 .48 .3 .1];
if nargin>=3
    h2.UserData.vol=vol;
    h2.UserData.opt=[];
elseif nargin<3
    btn_3D.Enable='inactive'; btn_3D.BackgroundColor=[1 1 1]*.9;
    btn_3D.ForegroundColor=[1 1 1]*.6;
    btn_3D.String='No 3D plot ';
end

%% Axes 4: Axial Slices - Upper Right Panel
% hp2 = uipanel(h,'Title','Axial','FontSize',12,...
%     'BackgroundColor','w','Foregroundcolor','k',...
%     'Position',[.51 .51 .48 .48],'Units','normalize');
ax(1)=axes('Position',[.51 .51 .45 .45]); box off;
imagesc(ax(1),squeeze(mri.anatomy(:,:,h2.UserData.z_axial_slice))); axis tight; view(-90,90); axis off;
hold on; p1=plot([0 mri.dim(2)],[h2.UserData.x_sagittal_slice h2.UserData.x_sagittal_slice],'r');
p2=plot([h2.UserData.y_coronal_slice h2.UserData.y_coronal_slice],[0 mri.dim(1)],'r');
vidx=h2.UserData.source_locs(:,3)==h2.UserData.z_axial_slice;
for v=fliplr(1:3)
    if vidx(v)==1
        a4(v)=scatter(h2.UserData.source_locs(v,2),h2.UserData.source_locs(v,1),'o','SizeData',h2.UserData.source_size,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    else    % plot off the image
        a4(v)=scatter(0,0,'o','SizeData',1,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    end
end
axis([0 h2.UserData.dims(2) 0 h2.UserData.dims(1)]);
%% Axes 3:  Sagittal Slices - Lower Left Panel
% hp3 = uipanel(h,'Title','Sagittal','FontSize',12,...
%     'BackgroundColor','w','Foregroundcolor','k',...
%     'Position',[.01 .01 .48 .48],'Units','normalize');
ax(2)=axes('Position',[.03 .03 .45 .45]); box off;
imagesc(ax(2),squeeze(mri.anatomy(h2.UserData.x_sagittal_slice,:,:))); axis tight; view(-90,90); axis off;
hold on; p1=plot([0 mri.dim(3)],[h2.UserData.y_coronal_slice h2.UserData.y_coronal_slice],'r');
p2=plot([h2.UserData.z_axial_slice h2.UserData.z_axial_slice],[0 mri.dim(2)],'r');
vidx=h2.UserData.source_locs(:,1)==h2.UserData.x_sagittal_slice;
for v=fliplr(1:3)
    if vidx(v)==1
        a3(v)=scatter(h2.UserData.source_locs(v,3),h2.UserData.source_locs(v,2),'o','SizeData',h2.UserData.source_size,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    else    % plot off the image
        a3(v)=scatter(0,0,'o','SizeData',1,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    end
end
axis([0 h2.UserData.dims(3) 0 h2.UserData.dims(2)]);
%% Axes 2:  Coronal Slices - Lower Right Panel
% hp4 = uipanel(h,'Title','Coronal','FontSize',12,...
%     'BackgroundColor','w','Foregroundcolor','k',...
%     'Position',[.51 .01 .48 .48],'Units','normalize');
ax(3)=axes('Position',[.51 .03 .45 .45]); box off;
imagesc(ax(3),squeeze(mri.anatomy(:,h2.UserData.y_coronal_slice,:))); axis tight; view(-90,90); axis off;
hold on; p1=plot([0 mri.dim(3)],[h2.UserData.x_sagittal_slice h2.UserData.x_sagittal_slice],'r');
p2=plot([h2.UserData.z_axial_slice h2.UserData.z_axial_slice],[0 mri.dim(1)],'r');
vidx=h2.UserData.source_locs(:,2)==h2.UserData.y_coronal_slice;
for v=fliplr(1:3)
    if vidx(v)==1
        a2(v)=scatter(h2.UserData.source_locs(v,3),h2.UserData.source_locs(v,1),'o','SizeData',h2.UserData.source_size,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    else    % plot off the image
        a2(v)=scatter(0,0,'o','SizeData',1,'MarkerFaceColor',h2.UserData.source_clr(v,:),'MarkerEdgeColor',h2.UserData.source_clr(v,:)*.3);
    end
end
axis([0 h2.UserData.dims(1) 0 h2.UserData.dims(3)]);

colormap(bone);

h2.Children(4).Children(end).Tag='axial'; h2.Children(4).Children(end).ButtonDownFcn=@startDragFcn_axial;
h2.Children(3).Children(end).Tag='sagittal'; h2.Children(3).Children(end).ButtonDownFcn=@startDragFcn_sagittal;
h2.Children(2).Children(end).Tag='coronal'; h2.Children(2).Children(end).ButtonDownFcn=@startDragFcn_coronal;
% for m=2:4; h2.Children(m).Children(3).PickableParts='all'; end

%% making source locs click & drag
for v=1:3
    h2.Children(4).Children(v).Tag=sprintf('axial Source %.f',v);
    h2.Children(4).Children(v).ButtonDownFcn=@startDragFcn_source;
    h2.Children(3).Children(v).Tag=sprintf('sagittal Source %.f',v);
    h2.Children(3).Children(v).ButtonDownFcn=@startDragFcn_source;
    h2.Children(2).Children(v).Tag=sprintf('coronal Source %.f',v);
    h2.Children(2).Children(v).ButtonDownFcn=@startDragFcn_source;
end

%% set up btn Down/Up funcs
h2.WindowButtonUpFcn=@stopDragFcn; % start ginput to get cursor location on axes
h2.WindowButtonDownFcn=@startDragFcn; % start ginput to get cursor location on axes

%% update image data for selected slice
function update_slices
global h2
if h2.UserData.x_sagittal_slice>0 && h2.UserData.x_sagittal_slice<=size(h2.UserData.anatomy,1) && ...
        h2.UserData.y_coronal_slice>0 && h2.UserData.y_coronal_slice<=size(h2.UserData.anatomy,2) && ...
        h2.UserData.z_axial_slice>0 && h2.UserData.z_axial_slice<=size(h2.UserData.anatomy,3)
    
    % update loc_text
    xyz=[h2.UserData.x_sagittal_slice,h2.UserData.y_coronal_slice,h2.UserData.z_axial_slice];
    xyz_w=ft_warp_apply(h2.UserData.transform,xyz); h2.UserData.current_loc=xyz_w;
    h2.UserData.source_locs_mm=ft_warp_apply(h2.UserData.transform,h2.UserData.source_locs); h2.UserData.current_loc=xyz_w;
    h2.UserData.loc_str=sprintf('pixel =   %.f   %.f   %.f\n mm =   %.f   %.f   %.f\n\nSource 1 =    %.f   %.f   %.f\nSource 2 =    %.f   %.f   %.f\nSource 3 =    %.f   %.f   %.f\n\nSource 1 =    %.f   %.f   %.f\nSource 2 =    %.f   %.f   %.f\nSource 3 =    %.f   %.f   %.f',...
        xyz,xyz_w,h2.UserData.source_locs',h2.UserData.source_locs_mm');
    h2.Children(1).Children(end-1).String=h2.UserData.loc_str; loc_txt.HorizontalAlignment='left';
    
    % updating images
    h2.Children(4).Children(end).CData=squeeze(h2.UserData.anatomy(:,:,h2.UserData.z_axial_slice)); % updating axial slice
    h2.Children(2).Children(end).CData=squeeze(h2.UserData.anatomy(:,h2.UserData.y_coronal_slice,:)); % updating axial slice
    h2.Children(3).Children(end).CData=squeeze(h2.UserData.anatomy(h2.UserData.x_sagittal_slice,:,:)); % updating axial slice
    
    % updating cursor lines
    % axial slice
    h2.Children(4).Children(end-2).XData=[0 size(h2.UserData.anatomy,2)];
    h2.Children(4).Children(end-2).YData=[h2.UserData.x_sagittal_slice h2.UserData.x_sagittal_slice];
    h2.Children(4).Children(end-1).XData=[h2.UserData.y_coronal_slice h2.UserData.y_coronal_slice];
    h2.Children(4).Children(end-1).YData=[0 size(h2.UserData.anatomy,1)];
    % sagital slice
    h2.Children(3).Children(end-2).XData=[h2.UserData.z_axial_slice h2.UserData.z_axial_slice];
    h2.Children(3).Children(end-2).YData=[0 size(h2.UserData.anatomy,2)];
    h2.Children(3).Children(end-1).YData=[h2.UserData.y_coronal_slice h2.UserData.y_coronal_slice];
    h2.Children(3).Children(end-1).XData=[0 size(h2.UserData.anatomy,3)];
    % coronal slice
    h2.Children(2).Children(end-2).XData=[h2.UserData.z_axial_slice h2.UserData.z_axial_slice];
    h2.Children(2).Children(end-2).YData=[0 size(h2.UserData.anatomy,1)];
    h2.Children(2).Children(end-1).YData=[h2.UserData.x_sagittal_slice h2.UserData.x_sagittal_slice];
    h2.Children(2).Children(end-1).XData=[0 size(h2.UserData.anatomy,3)];

    % updating Source locations when slice is selected
    % find all on axial slice
    vidx=h2.UserData.source_locs(:,3)>=h2.UserData.z_axial_slice-2 & h2.UserData.source_locs(:,3)<=h2.UserData.z_axial_slice+2; 
     for v=1:3
        if vidx(v)==1
            h2.Children(4).Children(v).XData=[h2.UserData.source_locs(v,2) h2.UserData.source_locs(v,2)];
            h2.Children(4).Children(v).YData=[h2.UserData.source_locs(v,1) h2.UserData.source_locs(v,1)];
            h2.Children(4).Children(v).SizeData=h2.UserData.source_size;
        else
           h2.Children(4).Children(v).XData=[0 0];
            h2.Children(4).Children(v).YData=[0 0];
            h2.Children(4).Children(v).SizeData=1;
        end
    end
    
    % find all on sagittal slice
    vidx=h2.UserData.source_locs(:,1)>=h2.UserData.x_sagittal_slice-2 & h2.UserData.source_locs(:,1)<=h2.UserData.x_sagittal_slice+2;
  for v=1:3
        if vidx(v)==1
            h2.Children(3).Children(v).XData=[h2.UserData.source_locs(v,3) h2.UserData.source_locs(v,3)];
            h2.Children(3).Children(v).YData=[h2.UserData.source_locs(v,2) h2.UserData.source_locs(v,2)];
            h2.Children(3).Children(v).SizeData=h2.UserData.source_size;
        else
           h2.Children(3).Children(v).XData=[0 0];
            h2.Children(3).Children(v).YData=[0 0];
            h2.Children(3).Children(v).SizeData=1;
        end
    end
    % find all on coronal slice
   vidx=h2.UserData.source_locs(:,2)>=h2.UserData.y_coronal_slice-2 & h2.UserData.source_locs(:,2)<=h2.UserData.y_coronal_slice+2;
   for v=1:3
        if vidx(v)==1
            h2.Children(2).Children(v).XData=[h2.UserData.source_locs(v,3) h2.UserData.source_locs(v,3)];
            h2.Children(2).Children(v).YData=[h2.UserData.source_locs(v,1) h2.UserData.source_locs(v,1)];
            h2.Children(2).Children(v).SizeData=h2.UserData.source_size;
        else
           h2.Children(2).Children(v).XData=[0 0];
            h2.Children(2).Children(v).YData=[0 0];
            h2.Children(2).Children(v).SizeData=1;
        end
    end
    
    
    
else
    h2.UserData.x_sagittal_slice=round(size(h2.UserData.anatomy,1)/2);
    h2.UserData.y_coronal_slice=round(size(h2.UserData.anatomy,2)/2);
    h2.UserData.z_axial_slice=round(size(h2.UserData.anatomy,3)/2);
end
function startDragFcn(varargin)
global h2
if strcmpi(varargin{2}.Source.CurrentAxes.Children(end).Tag,'axial')
    xyz=get(gca,'CurrentPoint');
    h2.UserData.x_sagittal_slice=round(xyz(3));
    h2.UserData.y_coronal_slice=round(xyz(1));
    update_slices;
    startDragFcn_axial;
elseif strcmpi(varargin{2}.Source.CurrentAxes.Children(end).Tag,'sagittal')
    xyz=get(gca,'CurrentPoint');
    h2.UserData.z_axial_slice=round(xyz(1));
    h2.UserData.y_coronal_slice=round(xyz(3));
    update_slices;
    startDragFcn_sagittal;
elseif strcmpi(varargin{2}.Source.CurrentAxes.Children(end).Tag,'coronal')
    xyz=get(gca,'CurrentPoint');
    h2.UserData.z_axial_slice=round(xyz(1));
    h2.UserData.x_sagittal_slice=round(xyz(3));
    update_slices;
    startDragFcn_coronal;
end
function stopDragFcn(varargin)
global h2
h2.WindowButtonMotionFcn='';   % stops dragging func
h2.UserData.drag_axial_source=0; % resetting so that sources are not selected to move during dragging
h2.UserData.drag_sagittal_source=0;
h2.UserData.drag_coronal_source=0;
function startDragFcn_source(varargin)
global h2
if strfind(varargin{2}.Source.Tag,'axial Source')==1
    h2.UserData.drag_axial_source=str2double(varargin{2}.Source.Tag(end));
    startDragFcn_axial;
elseif strfind(varargin{2}.Source.Tag,'sagittal Source')==1
    h2.UserData.drag_sagittal_source=str2double(varargin{2}.Source.Tag(end));
    startDragFcn_sagittal;
elseif strfind(varargin{2}.Source.Tag,'coronal Source')==1
   h2.UserData.drag_coronal_source=str2double(varargin{2}.Source.Tag(end));
    startDragFcn_coronal;
end

%% %%%%%%%%  Slice dragging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Axial Slice Dragging
function startDragFcn_axial(varargin)
global h2
h2.WindowButtonMotionFcn=@DragFcn_axial; % start ginput to get cursor location on axes
function DragFcn_axial(varargin)
global h2
xyz=get(gca,'CurrentPoint');
h2.UserData.x_sagittal_slice=round(xyz(3));
h2.UserData.y_coronal_slice=round(xyz(1));
if h2.UserData.drag_axial_source~=0
    v=h2.UserData.drag_axial_source;
    h2.UserData.source_locs(v,1)=round(xyz(3));
    h2.UserData.source_locs(v,2)=round(xyz(1));
end
update_slices;
%% Sagittal Slice Dragging
function startDragFcn_sagittal(varargin)
global h2
h2.WindowButtonMotionFcn=@DragFcn_sagittal; % start ginput to get cursor location on axes
function DragFcn_sagittal(varargin)
global h2
xyz=get(gca,'CurrentPoint');
h2.UserData.z_axial_slice=round(xyz(1));
h2.UserData.y_coronal_slice=round(xyz(3));
if h2.UserData.drag_sagittal_source~=0
    v=h2.UserData.drag_sagittal_source;
    h2.UserData.source_locs(v,3)=round(xyz(1));
    h2.UserData.source_locs(v,2)=round(xyz(3));
end
update_slices;
%% Coronal Slice Dragging
function startDragFcn_coronal(varargin)
global h2
h2.WindowButtonMotionFcn=@DragFcn_coronal; % start ginput to get cursor location on axes
function DragFcn_coronal(varargin)
global h2
xyz=get(gca,'CurrentPoint');
h2.UserData.z_axial_slice=round(xyz(1));
h2.UserData.x_sagittal_slice=round(xyz(3));
if h2.UserData.drag_coronal_source~=0
    v=h2.UserData.drag_coronal_source;
    h2.UserData.source_locs(v,3)=round(xyz(1));
    h2.UserData.source_locs(v,1)=round(xyz(3));
end
update_slices;


%% Select Source location
function select_location(varargin)
global h2
if strcmpi(varargin{2}.Source.String,'Set Source 1')
    h2.UserData.source_locs(1,:)=[h2.UserData.x_sagittal_slice h2.UserData.y_coronal_slice h2.UserData.z_axial_slice];
elseif strcmpi(varargin{2}.Source.String,'Set Source 2')
    h2.UserData.source_locs(2,:)=[h2.UserData.x_sagittal_slice h2.UserData.y_coronal_slice h2.UserData.z_axial_slice];
elseif strcmpi(varargin{2}.Source.String,'Set Source 3')
    h2.UserData.source_locs(3,:)=[h2.UserData.x_sagittal_slice h2.UserData.y_coronal_slice h2.UserData.z_axial_slice];
end
update_slices;
function goto_location(varargin)
global h2
if strcmpi(varargin{2}.Source.String,'Go to Source 1')
    h2.UserData.x_sagittal_slice=h2.UserData.source_locs(1,1);
    h2.UserData.y_coronal_slice=h2.UserData.source_locs(1,2);
    h2.UserData.z_axial_slice=h2.UserData.source_locs(1,3);
elseif strcmpi(varargin{2}.Source.String,'Go to Source 2')
    h2.UserData.x_sagittal_slice=h2.UserData.source_locs(2,1);
    h2.UserData.y_coronal_slice=h2.UserData.source_locs(2,2);
    h2.UserData.z_axial_slice=h2.UserData.source_locs(2,3);
elseif strcmpi(varargin{2}.Source.String,'Go to Source 3')
    h2.UserData.x_sagittal_slice=h2.UserData.source_locs(3,1);
    h2.UserData.y_coronal_slice=h2.UserData.source_locs(3,2);
    h2.UserData.z_axial_slice=h2.UserData.source_locs(3,3);
end
update_slices;
%% save source locs to workspace
function save_location(varargin)
global h2 h
source_locs=h2.UserData.source_locs;
source_locs_mm=h2.UserData.source_locs_mm;
% assignin('base','source_locs',source_locs);
% assignin('base','source_locs_mm',source_locs_mm);
h.cfg.study.source_locs=h2.UserData.source_locs;
h.cfg.study.source_locs_mm=h2.UserData.source_locs_mm;

%% Plot source_locs_mm in 3D volume "vol"
function plot_source_locs_3D(varargin)
global h2
if ~isfield(h2.UserData.opt,'axes_h')
    figure; set(gcf,'color','w'); cla; ax=gca; axis off;
    h2.UserData.opt.axes_h=ax;
else
    if ~isvalid(h2.UserData.opt.axes_h)
        figure; set(gcf,'color','w'); cla; ax=gca;  axis off;
        h2.UserData.opt.axes_h=ax;
    end
end
h2.UserData.opt.source_clr=h2.UserData.source_clr;
h2.UserData.opt.vol_nums=1:length(h2.UserData.vol);
bl_plot_source_locs(h2.UserData.vol,h2.UserData.source_locs_mm,h2.UserData.opt); 

function default_location(varargin)
global h2
h2.UserData.source_locs=[35 109 78; 151 109 78; 91 48 78]; %ones(3,3); %[];
update_slices();
