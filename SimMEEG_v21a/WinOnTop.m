function wasOnTop = WinOnTop( figureHandle, isOnTop )
%WINONTOP allows to trigger figure's "Always On Top" state
%
%% INPUT ARGUMENTS:
%
% # figureHandle - Matlab's figure handle, scalar
% # isOnTop      - logical scalar or empty array
%
%
%% USAGE:
%
% * WinOnTop( hfigure, true );      - switch on  "always on top"
% * WinOnTop( hfigure, false );     - switch off "always on top"
% * WinOnTop( hfigure );            - equal to WinOnTop( hfigure,true);
% * WinOnTop();                     - equal to WinOnTop( gcf, true);
% * WasOnTop = WinOnTop(...);       - returns boolean value "if figure WAS on top"
% * isOnTop = WinOnTop(hfigure,[])  - get "if figure is on top" property
%
% For Matlab windows, created via `hf=uifigure()` use `uifigureOnTop()`, see: 
% https://www.mathworks.com/matlabcentral/fileexchange/73134-uifigureontop
%
%% LIMITATIONS:
%
% * java enabled
% * figure must be visible
% * figure's "WindowStyle" should be "normal"
% * figureHandle should not be casted to double, if using HG2 (R2014b+)
%
%
% Written by Igor
% i3v@mail.ru
%
% 2013.06.16 - Initial version
% 2013.06.27 - removed custom "ishandle_scalar" function call
% 2015.04.17 - adapted for changes in matlab graphics system (since R2014b)
% 2016.05.21 - another ishg2() checking mechanism 
% 2016.09.24 - fixed IsOnTop vs isOnTop bug
% 2019.10.27 - link for uifigureOnTop; connected to github; renamed to "demo_"
% 
% Copyright (c) 2016, Igor
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.



%% Parse Inputs

if ~exist('figureHandle','var'); figureHandle = gcf; end

assert(...
          isscalar(  figureHandle ) &&...
          ishandle(  figureHandle ) &&...
          strcmp(get(figureHandle,'Type'),'figure'),...
          ...
          'WinOnTop:Bad_figureHandle_input',...
          '%s','Provided figureHandle input is not a figure handle'...
       );

assert(...
            strcmp('on',get(figureHandle,'Visible')),...
            'WinOnTop:FigInisible',...
            '%s','Figure Must be Visible'...
       );

assert(...
            strcmp('normal',get(figureHandle,'WindowStyle')),...
            'WinOnTop:FigWrongWindowStyle',...
            '%s','WindowStyle Must be Normal'...
       );
   
if ~exist('isOnTop','var'); isOnTop=true; end

assert(...
          islogical( isOnTop ) && ...
          isscalar(  isOnTop ) || ...
          isempty(   isOnTop ),  ...
          ...
          'WinOnTop:Bad_isOnTop_input',...
          '%s','Provided isOnTop input is neither boolean, nor empty'...
      );
  
  
%% Pre-checks

error(javachk('swing',mfilename)) % Swing components must be available.
  
  
%% Action

% Flush the Event Queue of Graphic Objects and Update the Figure Window.
drawnow expose

warnStruct=warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(handle(figureHandle),'JavaFrame');
warning(warnStruct.state,'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

drawnow


if ishg2(figureHandle)
    jFrame_fHGxClient = jFrame.fHG2Client;
else
    jFrame_fHGxClient = jFrame.fHG1Client;
end


wasOnTop = jFrame_fHGxClient.getWindow.isAlwaysOnTop;

if ~isempty(isOnTop)
    jFrame_fHGxClient.getWindow.setAlwaysOnTop(isOnTop);
end

end


function tf = ishg2(figureHandle)
% There's a detailed discussion, how to check "if using HG2" here:
% http://www.mathworks.com/matlabcentral/answers/136834-determine-if-using-hg2
% however, it looks like there's no perfect solution.
%
% This approach, suggested by Cris Luengo:
% http://www.mathworks.com/matlabcentral/answers/136834#answer_156739
% should work OK, assuming user is NOT passing a figure handle, casted to
% double, like this:
%
%   hf=figure();
%   WinOnTop(double(hf));
%

tf = isa(figureHandle,'matlab.ui.Figure');

end