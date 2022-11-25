function d = dim(Z)
% dim - returns the dimension in which the zonotope is defined;
%    caution: this is different from the rank of a zonotope
%
% Syntax:  
%    d = dim(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    d - dimension in which Z is defined
%
% Example: 
%    Z = zonotope([zeros(3,1),rand(3,5)]);
%    d = dim(Z)
%
% Other m-files required: center.m
% Subfunctions: none
% MAT-files required: none
%
% See also: rank.m
%
% Author:        Mark Wetzlinger
% Written:       15-Sep-2019
% Last update:   11-March-2021 (MW, add empty case)
% Last revision: ---

%------------- BEGIN CODE --------------

if ~isempty(Z)
    d = length(center(Z));
else
    d = 0;
end

%------------- END OF CODE --------------