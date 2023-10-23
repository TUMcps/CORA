function c = center(E)
% center - returns the center of an ellipsoid object
%
% Syntax:
%    c = center(E);
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    c - center of ellipsoid
%
% Example: 
%    E = ellipsoid([1 0; 0 1]);
%    c = center(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   04-July-2022 (VG, avoid class array problems)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% make sure it is not a class array
inputArgsCheck({{E,'att','ellipsoid','scalar'}});

c = E.q;

% ------------------------------ END OF CODE ------------------------------
