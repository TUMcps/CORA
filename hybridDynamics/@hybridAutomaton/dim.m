function n = dim(HA)
% dim - returns the dimension of each location; if a single hybrid
%    automaton is given and all locations have the same dimension, only one
%    value is returned
%
% Syntax:  
%    n = dim(HA)
%
% Inputs:
%    HA - hybridAutomaton object
%
% Outputs:
%    n - vector of dimensions
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      22-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check size
[r,c] = size(HA);

if r > 1 && c > 1
    throw(CORAerror('CORA:notSupported',...
        'Only supported for a row/column array of hybridAutomaton objects.'));
end

% number of hybrid automata
numHA = length(HA);

if numHA > 1
    % use cell-array since HA may have different number of locations
    n = cell(numHA,1);
    for i=1:numHA
        n{i} = arrayfun(@(x) x.contDynamics.dim,HA(i).location,...
            'UniformOutput',true);
    end

else

    n = arrayfun(@(x) x.contDynamics.dim,HA.location,'UniformOutput',true);
    if all(n == n(1))
        n = n(1);
    end

end

%------------- END OF CODE --------------