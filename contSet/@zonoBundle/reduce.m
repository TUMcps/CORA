function zB = reduce(zB,option,varargin)
% reduce - Reduces the order of a zonotope bundle
%
% Syntax:
%    zB = reduce(zB,option,order,filterLength)
%
% Inputs:
%    zB - zonotope bundle
%    option - reduction method selector
%    order - maximum order of reduced zonotope
%    filterLength - ???
%
% Outputs:
%    zB - bundle of reduced zonotopes
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce

% Authors:       Matthias Althoff
% Written:       09-November-2010 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% reduce order of each zonotope
for i=1:zB.parallelSets
    zB.Z{i} = reduce(zB.Z{i},option,varargin{:});
end

% ------------------------------ END OF CODE ------------------------------
