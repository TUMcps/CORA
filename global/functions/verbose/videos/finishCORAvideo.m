function res = finishCORAvideo(vidObj,fig)
% finishCORAvideo - finishes the CORA video
%
% Syntax:
%    res = finishCORAvideo(vidObj,fig)
%
% Inputs:
%    vidObj - VideoWriter
%    fig - figure
%
% Outputs:
%    res - logical
%
% See also:
%    CORAvideo_snippets

% Authors:       Tobias Ladner
% Written:       05-April-2024
% Last update:   11-December-2024 (TL, split setup and recording)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp('Finishing up CORA video..')

% close everything
close(vidObj)
close(fig)

videoName = [vidObj.Path filesep vidObj.Filename];
fprintf('- Video saved to: <a href="matlab:winopen(''%s'');">%s</a> \n', vidObj.Path, videoName)

% completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
