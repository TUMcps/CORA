function [val,x] = supportFunc(obj,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a polytope along a
%    certain direction
%
% Syntax:  
%    val = supportFunc(obj,dir)
%    [val,x] = supportFunc(obj,dir,type)
%
% Inputs:
%    obj - mptPolytope object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the polytope in the specified direction
%    x - support vector
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/supportFunc

% Author:       Niklas Kochdumper, Victor Gassmann
% Written:      19-November-2019
% Last update:  16-March-2021 (added unbounded support)
% Last revision:---

%------------- BEGIN CODE --------------
    
% pre-processing
[res,vars] = pre_supportFunc('mptPolytope',obj,dir,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    val = vars{1}; return
else
    obj = vars{1}; dir = vars{2}; type = vars{3};
end


% upper or lower bound
if strcmp(type,'lower')
    s = 1;
elseif strcmp(type,'upper')
    s = -1;
elseif strcmp(type,'range')
    throw(CORAerror('CORA:notSupported',type));
end

% get object properties
A = obj.P.A;
b = obj.P.b;

% linear program options
options = optimoptions('linprog','display','off');

[x,val,exitflag] = linprog(s*dir',A,b,[],[],[],[],options);
val = s*val;

if exitflag == -3
    % unbounded
    val = -s*Inf;
    x = -s*sign(dir).*Inf(length(dir),1);
elseif exitflag ~= 1
    throw(CORAerror('CORA:solverIssue'));
end

%------------- END OF CODE --------------