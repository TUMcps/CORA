function val = norm(E,varargin)
% norm - compute the maximum Euclidean norm of an ellipsoid
%
% Syntax:
%    val = norm(E)
%    val = norm(E,type)
%
% Inputs:
%    E    - ellipsoid object 
%    type - (optional) norm type (default: 2)
%
% Outputs:
%    val - value of the maximum norm
%
% Example: 
%    E = ellipsoid.generateRandom('Dimension',2,'Center',zeros(2,1));
%    val = norm(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      20-November-2019
% Last update:  31-July-2020
%               04-July-2022 (VG: class array case)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,var] = pre_norm('ellipsoid',E,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    val = var{1}; return
else
    % assign values
    E = var{1}; type = var{2};
end


% only Euclidean norm implemented
if ~isnumeric(type) || type ~= 2
    throw(CORAerror('CORA:noExactAlg',E,type,...
        'Only implemented for Euclidean norm'));
end

% only for zero-centers implemented
if ~all(E.q==0)
    throw(CORAerror('CORA:noExactAlg',E,type,...
        'Not yet implemented for non-zero center'));
end

% transform into eigenspace
lmax = max(eig(E.Q));
val = sqrt(lmax);

%------------- END OF CODE --------------