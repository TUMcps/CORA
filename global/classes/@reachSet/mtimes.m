function R = mtimes(factor1,factor2)
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
% See also: contSet/mtimes

% Authors:       Niklas Kochdumper
% Written:       04-March-2021
% Last update:   10-November-2022 (MW, add checks for empty structs)
%                02-February-2024 (TL, use findClassArgs for R * 2 to work)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

[R,M] = findClassArg(factor1,factor2,'reachSet');

for i = 1:size(R,1)
    if ~isempty(R(i).timeInterval)
        R(i).timeInterval.set = cellfun(@(S) M*S, ...
           R(i).timeInterval.set,'UniformOutput',false);
    end
    if ~isempty(R(i).timePoint)
        R(i).timePoint.set = cellfun(@(S) M*S, ...
            R(i).timePoint.set,'UniformOutput',false);
    end
end

% ------------------------------ END OF CODE ------------------------------
