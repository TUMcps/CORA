function R = mtimes(M,R)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%    with a reachSet object
%
% Syntax:
%    R = mtimes(M,R)
%
% Inputs:
%    M - numerical matrix
%    R - reachSet object 
%
% Outputs:
%    R - transformed reachSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Niklas Kochdumper
% Written:       04-March-2021
% Last update:   10-November-2022 (MW, add checks for empty structs)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for i = 1:size(R,1)
    if ~isempty(R(i).timeInterval)
        R(i).timeInterval.set = cellfun(@(x) M*x, ...
           R(i).timeInterval.set,'UniformOutput',false);
    end
    if ~isempty(R(i).timePoint)
        R(i).timePoint.set = cellfun(@(x) M*x, ...
            R(i).timePoint.set,'UniformOutput',false);
    end
end

% ------------------------------ END OF CODE ------------------------------
