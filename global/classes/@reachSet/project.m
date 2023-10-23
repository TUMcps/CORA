function R = project(R,dims)
% project - projects reachable set in reachSet object to specified
%    dimensions
%
% Syntax:
%    R = project(R,dims)
%
% Inputs:
%    R - reachSet object 
%    dims - dimensions onto which reachSet should be projected
%
% Outputs:
%    R - resulting transformed reachset object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Mark Wetzlinger
% Written:       18-November-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% project 
for i = 1:size(R,1)

    % time-point reachable set
    if ~isempty(R(i).timePoint)
        R(i).timePoint.set = cellfun(@(x) project(x,dims), ...
            R(i).timePoint.set,'UniformOutput',false);
    end

    % time-interval reachable set
    if ~isempty(R(i).timeInterval)
        R(i).timeInterval.set = cellfun(@(x) project(x,dims), ...
            R(i).timeInterval.set,'UniformOutput',false);
    end
end

% ------------------------------ END OF CODE ------------------------------
