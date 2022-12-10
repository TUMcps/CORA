function [val,x] = supportFunc(C,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a capsule along a
%    certain direction
%
% Syntax:  
%    [val,x] = supportFunc(C,dir)
%    [val,x] = supportFunc(C,dir,type)
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
% See also: conZonotope/supportFunc

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  10-December-2022 (MW, add type = 'range')
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_supportFunc('capsule',C,dir,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    val = vars{1}; x = []; return
else
    C = vars{1}; dir = vars{2}; type = vars{3};
end


% get object properties
c = C.c;
g = C.g;
r = C.r;

% consider case where capsule is just a ball
if isempty(g)
    g = zeros(size(c)); 
end

% compute upper or lower bound
if strcmp(type,'upper')
    val = dir'*c + abs(dir'*g) + r*norm(dir);  
    x = c + g*sign(dir'*g) + r*dir;
elseif strcmp(type,'lower')
    val = dir'*c - abs(dir'*g) - r*norm(dir);  
    x = c - g*sign(dir'*g) - r*dir;
elseif strcmp(type,'range')
    val = interval(dir'*c - abs(dir'*g) - r*norm(dir),...
        dir'*c + abs(dir'*g) + r*norm(dir));
    x = [c - g*sign(dir'*g) - r*dir, c + g*sign(dir'*g) + r*dir];
end

%------------- END OF CODE --------------