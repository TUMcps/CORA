function val = norm(obj,type)
% norm - computes the exact maximum norm value of specified norm
%
% Syntax:  
%    val = norm(obj,varargin)
%
% Inputs:
%    obj - interval object
%    type - additional arguments of builtin/norm
%
% Outputs:
%    val - norm value
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Victor Gassmann
% Written:      31-July-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('type','var')
    type = 2;
end
if  ~isnumeric(type) && ~strcmp(type,'fro')
    error('Wrong type of argument for norm type');
end
if isnumeric(type) && type==2
    val = norm(max(abs(obj.inf),abs(obj.sup)));
else
    val = norm(intervalMatrix(center(obj),rad(obj)),type);
end
%------------- END OF CODE --------------