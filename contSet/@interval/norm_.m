function val = norm_(I,type,varargin)
% norm_ - computes the exact maximum norm value of specified norm
%
% Syntax:
%    val = norm_(I,type)
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
% See also: contSet/norm, mtimes

% Authors:       Victor Gassmann
% Written:       31-July-2020
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename norm_)

% ------------------------------ BEGIN CODE -------------------------------

if isnumeric(type) && type == 2
    if isempty(I.inf)
        % this exception handling for empty intervals is required, since
        % MATLAB computes norm([]) = 0, but we need -Inf
        val = -Inf;
    else
        val = norm(max(abs(I.inf),abs(I.sup)));
    end
else
    val = norm(intervalMatrix(center(I),rad(I)),type);
end

% ------------------------------ END OF CODE ------------------------------
