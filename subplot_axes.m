function [ax]=subplot_axes(num_rows,num_col,x_sep,y_sep,x_overlap,y_overlap,plot_order)
%function [ax]=subplot_axes(num_rows,num_col,x_sep,y_sep,x_overlap,y_overlap,plot_order)
% INPUT
%   num_rows = number fo rows
%   num_col = number of columns
%   x_sep = horizontal separation between plots considering x-scale = [0 1].
%   y_sep = vertical separation between plots considering y-scale = [0 1].
%   x_overlap = horizontal overlap between plots considering, must be < x_sep.
%   y_overlap = vertical overlap between plots considering, must be < y_sep.
%   plot_order = plot in order along (0) columns then rows [default] or (1) rows then columns
%
% OUTPUT
%   ax = index locator for all axes created. starts in bottom left corner going left-to-right then up by rows.
%    Note: axis off by default, bu text labelled for convenience
%
% Created by A. Herdman to replace matlab's subplot.m function.
% Feb 8, 2015
if nargin<5; x_overlap=0; y_overlap=0; plot_order=0;
elseif nargin<6; y_overlap=0; plot_order=0;
elseif nargin<7; plot_order=0;
end


x_size=(1-x_sep*(num_col+1))/num_col;
y_size=(1-y_sep*(num_rows+1))/num_rows;

x = [x_sep:x_size+x_sep-x_overlap:1];%+((num_rows-1)*-x_overlap)
y = [y_sep:y_size+y_sep-y_overlap:1]+((num_rows-1)*y_overlap);

x=x(1:num_col);
y=fliplr(y(1:num_rows));
a=1;
if plot_order==0 % default order ==> axes across columns then rows
    for v=1:num_rows
        for h=1:num_col
            ax(a)=axes('position',[x(h) y(v) x_size y_size]);
 %           text(0.5,0.5,num2str(a));
            axis off;
            a=a+1;
        end
    end
elseif plot_order==1 % default order ==> axes across columns then rows
    for h=1:num_col
        for v=1:num_rows
            ax(a)=axes('position',[x(h) y(v) x_size y_size]);
%             text(0.5,0.5,num2str(a));
            axis off;
            a=a+1;
        end
    end
end



