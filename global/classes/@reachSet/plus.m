function R = plus(R,S)
% plus - Overloaded '+' operator for the Minkowski addition of a set or a
%    vector with a reachSet object
%
% Syntax:
%    R = plus(R,S)
%
% Inputs:
%    R - reachSet object 
%    S - contSet object or numerical vector
%
% Outputs:
%    R - resulting tranformed reachset object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Niklas Kochdumper
% Written:       04-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get reachSet object
[R,S] = findClassArg(R,S,'reachSet');

% compute Minkowski sum
for i = 1:size(R,1)
    if ~isempty(R(i).timeInterval)
        R(i).timeInterval.set = cellfun(@(RS) RS + S,...
                        R(i).timeInterval.set,'UniformOutput',false);
    end
    if ~isempty(R(i).timePoint)
        R(i).timePoint.set = cellfun(@(RS) RS + S,...
                        R(i).timePoint.set,'UniformOutput',false); 
    end
end

% ------------------------------ END OF CODE ------------------------------
