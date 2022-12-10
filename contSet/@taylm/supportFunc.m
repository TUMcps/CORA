function val = supportFunc(tay,dir,varargin)
% supportFunc - Calculate the upper or lower bound of a Taylor model along 
%    a certain direction
%
% Syntax:  
%    val = supportFunc(tay,dir)
%    val = supportFunc(tay,dir,type)
%
% Inputs:
%    tay - taylm object
%    dir - direction for which the bounds are calculated (vector of size
%          (n,1) )
%    type - upper or lower bound ('lower' or 'upper')
%
% Outputs:
%    val - bound of the Taylor model in the specified direction
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/supportFunc

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
% pre-processing
[res,vars] = pre_supportFunc('taylm',tay,dir,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    val = vars{1}; return
else
    tay = vars{1}; dir = vars{2}; type = vars{3};
end

% compute enclosing interval of Taylor model projeted onto the direction
I = interval(dir'*tay);

% upper or lower bound
if strcmp(type,'lower')
    val = infimum(I);
elseif strcmp(type,'range')
    val = supremum(I);
else
    % 'val' is already the desired result
end

%------------- END OF CODE --------------