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
%    type - upper bound, lower bound, or both ('upper','lower','range')
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
%               10-December-2022 (MW, add 'range')
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


% linear program options
options = optimoptions('linprog','display','off');

% upper or lower bound
if strcmp(type,'lower')
    [val,x] = aux_solveLinProg(obj,dir,1,options);
elseif strcmp(type,'upper')
    [val,x] = aux_solveLinProg(obj,dir,-1,options);
elseif strcmp(type,'range')
    [val_upper,x_upper] = aux_solveLinProg(obj,dir,1,options);
    [val_lower,x_lower] = aux_solveLinProg(obj,dir,-1,options);
    % combine values
    val = interval(val_lower,val_upper);
    x = [x_lower x_upper];
end

end


% Auxiliary function ------------------------------------------------------

function [val,x] = aux_solveLinProg(obj,dir,s,options)

% solve linear program
[x,val,exitflag] = linprog(s*dir',obj.P.A,obj.P.b,[],[],[],[],options);
val = s*val;

if exitflag == -3
    % unbounded
    val = -s*Inf;
    x = -s*sign(dir).*Inf(length(dir),1);
elseif exitflag ~= 1
    throw(CORAerror('CORA:solverIssue'));
end

end

%------------- END OF CODE --------------