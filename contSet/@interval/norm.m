function val = norm(I,varargin)
% norm - computes the exact maximum norm value of specified norm
%
% Syntax:  
%    val = norm(I,varargin)
%
% Inputs:
%    I - interval object
%    type - (optional) additional arguments of builtin/norm
%
% Outputs:
%    val - norm value
%
% Example:
%    I = interval([-2;1],[3;2]);
%    norm(I,2)
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

% pre-processing
[res,var] = pre_norm('interval',I,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    val = var{1}; return
else
    % assign values
    I = var{1}; type = var{2};
end


if isnumeric(type) && type == 2
    val = norm(max(abs(I.inf),abs(I.sup)));
else
    val = norm(intervalMatrix(center(I),rad(I)),type);
end

%------------- END OF CODE --------------