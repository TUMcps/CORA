function R = times(v,R)
% times - Overloaded '.*' operator for the multiplication of a matrix or an
%    with a reachSet object
%
% Syntax:
%    R = times(v,R)
%
% Inputs:
%    v - numerical vector
%    R - reachSet object 
%
% Outputs:
%    R - transformed reachSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/times

% Authors:       Tobias Ladner
% Written:       02-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for i = 1:size(R,1)
    if ~isempty(R(i).timeInterval)
        R(i).timeInterval.set = cellfun(@(x) v.*x, ...
           R(i).timeInterval.set,'UniformOutput',false);
    end
    if ~isempty(R(i).timePoint)
        R(i).timePoint.set = cellfun(@(x) v.*x, ...
            R(i).timePoint.set,'UniformOutput',false);
    end
end

% ------------------------------ END OF CODE ------------------------------
