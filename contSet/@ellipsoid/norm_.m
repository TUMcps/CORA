function val = norm_(E,type,varargin)
% norm_ - compute the maximum Euclidean norm of an ellipsoid
%
% Syntax:
%    val = norm_(E)
%    val = norm_(E,type)
%
% Inputs:
%    E    - ellipsoid object 
%    type - (optional) norm type (default: 2)
%
% Outputs:
%    val - value of the maximum norm
%
% Example: 
%    E = ellipsoid([3 -1; -1 1],[0;0]);
%    val = norm(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/norm

% Authors:       Victor Gassmann
% Written:       20-November-2019
% Last update:   31-July-2020
%                04-July-2022 (VG, class array case)
% Last revision: 27-March-2023 (MW, rename norm_)

% ------------------------------ BEGIN CODE -------------------------------

% only Euclidean norm implemented
if ~isnumeric(type) || type ~= 2
    throw(CORAerror('CORA:noExactAlg',E,type,...
        'Only implemented for Euclidean norm'));
end

% only for zero-centers implemented
if ~all(E.q==0)
    throw(CORAerror('CORA:notSupported',...
        'Not yet implemented for non-zero center.'));
end

% transform into eigenspace
lmax = max(eig(E.Q));
val = sqrt(lmax);

% check for empty set
if isnan(val) && representsa(E,'emptySet',0)
    val = -Inf;
end

% ------------------------------ END OF CODE ------------------------------
