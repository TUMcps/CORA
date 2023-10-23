function [val,x] = supportFunc_(C,dir,type,varargin)
% supportFunc_ - Calculate the upper or lower bound of a capsule along a
%    certain direction
%
% Syntax:
%    [val,x] = supportFunc_(C,dir)
%    [val,x] = supportFunc_(C,dir,type)
%
% Inputs:
%    C - capsule object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper bound, lower bound, or both ('upper','lower','range')
%
% Outputs:
%    val - bound of the capsule in the specified direction (if type =
%          'upper' or 'lower'), or interval with both (type = 'range')
%    x - support vector(s)
%
% Example: 
%    C = capsule([1; 1], [0.5; -1], 0.5);
%    val = supportFunc(C,[1;1]);
%   
%    figure; hold on;
%    plot(C,[1,2],'r');
%    plot(conHyperplane(halfspace([1;1],val)),[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/supportFunc, conZonotope/supportFunc_

% Authors:       Niklas Kochdumper
% Written:       19-November-2019
% Last update:   10-December-2022 (MW, add type = 'range')
%                25-April-2023 (MW, bug fix)
% Last revision: 27-March-2023 (MW, rename supportFunc_)

% ------------------------------ BEGIN CODE -------------------------------

% get object properties
c = C.c;
g = C.g;
r = C.r;

% consider case where capsule is just a ball
if isempty(g)
    g = zeros(size(c)); 
end

% length of direction
l = vecnorm(dir);

% compute upper or lower bound
if strcmp(type,'upper')
    val = dir'*c + abs(dir'*g) + r*l;  
    x = c + g*sign(dir'*g) + r * dir ./ l;

elseif strcmp(type,'lower')
    val = dir'*c - abs(dir'*g) - r*l;  
    x = c - g*sign(dir'*g) - r * dir ./ l;

elseif strcmp(type,'range')
    val = interval(dir'*c - abs(dir'*g) - r*l,...
        dir'*c + abs(dir'*g) + r*l);
    x = [c - g*sign(dir'*g) - r*dir./l, c + g*sign(dir'*g) + r*dir./l];

end

% ------------------------------ END OF CODE ------------------------------
