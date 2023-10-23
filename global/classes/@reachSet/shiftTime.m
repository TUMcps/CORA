function R = shiftTime(R,delta)
% shiftTime - shifts all sets of a reachSet object by a scalar time delta
%
% Syntax:
%    R = shiftTime(R,delta)
%
% Inputs:
%    R - reachSet object 
%    delta - scalar
%
% Outputs:
%    R - shifted reachset object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Mark Wetzlinger
% Written:       24-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get reachSet object
[R,delta] = findClassArg(R,delta,'reachSet');

% shift time vector
for i = 1:size(R,1)
    if ~isempty(R(i).timeInterval)
        R(i).timeInterval.time = cellfun(@(x) delta + x,...
                        R(i).timeInterval.time,'UniformOutput',false);
    end
    if ~isempty(R(i).timePoint)
        R(i).timePoint.time = cellfun(@(x) delta + x,...
                        R(i).timePoint.time,'UniformOutput',false); 
    end
end

% ------------------------------ END OF CODE ------------------------------
