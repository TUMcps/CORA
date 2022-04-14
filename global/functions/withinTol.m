function res = withinTol(a,b,TOL)
% sumPoints - checks whether a,b are within a given tolerance
%
% Syntax:  
%    res = withinTol(a,b)
%    res = withinTol(a,b,TOL)
%
% Inputs:
%    a,b - double to check tolerances for
%
% Outputs:
%    res - boolean result
%
% Example: 
%    res = withinTol(1,1+1e-12)
%
% Other m-files required: -
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      19-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('TOL','var')
    TOL = 1e-8;
end

% allow scalar values to be expanded
if ~all(size(a)==size(b)) 
    if isscalar(a)
        a = repmat(a,size(b));
    elseif isscalar(b)
        b = repmat(b,size(a));
    else
        error('First two arguments need to be of same size (or scalar)!');
    end
end

if ~isa(a,'double') || ~isa(b,'double')
    error('First two input arguments need to be doubles!');
end

% absolute tolerance
res_abs = abs(a-b)<=TOL;

% relative tolerance
res_rel = abs(a-b)./min(abs(a),abs(b))<=TOL;

res = res_abs | res_rel;
%------------- END OF CODE --------------