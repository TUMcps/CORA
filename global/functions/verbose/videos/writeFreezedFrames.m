function writeFreezedFrames(vidObj,duration)
% writeFreezedFrames - writes frames such that the video stays for a few
%   seconds
%
% Syntax:
%    writeFreezedFrames(vidObj,duration)
%
% Inputs:
%    vidObj - VideoWriter
%    duration - numeric positive scalar
%
% Outputs:
%    -
%
% See also:
%    setupCORAvideo

% Authors:       Tobias Ladner
% Written:       10-December-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if duration > 0
    % read out frame rate
    FrameRate = vidObj.FrameRate;
    
    % write frames
    writeVideo(vidObj, repmat(getframe(gcf),FrameRate*duration,1));
end

end

% ------------------------------ END OF CODE ------------------------------
