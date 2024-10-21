function sig = cutAtFirstFallingEdge(obj)
% cutAtFirstFallingEdge - cut the signal at the first falling edge, i.e., the first time the signal changes from true to false
%                         all values after the falling edge are set to false
%
% Syntax:
%    sig = cutAtFirstFallingEdge(obj)
%
% Inputs:
%    obj - pointSegmentSignal object
%
% Outputs:
%    sig - the cut signal
%
% Example:
%    obj = pointSegmentSignal([0 2 3],[true true false false false true]);
%    sig = cutAtFirstFallingEdge(obj)
%    assert(sig == pointSegmentSignal([0 2],[true true false false]));
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       12-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% find the first time with a falling edge
i = find(~obj.values(2:end) & obj.values(1:end-1),1,'first');
if ~isempty(i)
    if mod(i,2) ~= 0
        % falls from point to interval
        sig = pointSegmentSignal(obj.timePoints(1:((i + 1) / 2)),[obj.values(1:i),false]);
    else
        % falls from interval to point
        sig = pointSegmentSignal(obj.timePoints(1:((i / 2) + 1)),[obj.values(1:i),false,false]);
    end
else
    % there is no falling edge, i.e., there is nothing to cut
    % thus, we return the original signal
    sig = obj;
end

% ------------------------------ END OF CODE ------------------------------
