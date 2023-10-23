function E = project(E,dims)
% project - projects an ellipsoid onto the specified dimensions
%
% Syntax:
%    E = project(E,dims)
%
% Inputs:
%    E - (ellipsoid) ellipsoid
%    dims - dimensions for projection
%
% Outputs:
%    E - (ellipsoid) projected ellipsoid
%
% Example: 
%    E = ellipsoid([9.3 -0.6 1.9;-0.6 4.7 2.5; 1.9 2.5 4.2]);
%    E = project(E,[1 3])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   04-July-2022 (VG, input checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments
inputArgsCheck({{E,'att','ellipsoid','scalar'};
                {dims,'att',{'numeric','logical'},{{'nonnan','vector',...
                     'integer','>=',1,'<=',dim(E)},{'vector'}}}});

% project set
I = eye(length(E.Q));
P = I(:,dims);
E = P'*E;

% ------------------------------ END OF CODE ------------------------------
