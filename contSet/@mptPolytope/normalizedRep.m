function P = normalizedRep(P)
% normalizedRep - Computes a normalized representation by scaling the
%                 length of the normal vectors to 1
%
% Syntax:  
%    P = normalizedRep(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    P - mptPolytope object
%
% Example: 
%    P = mptPolytope([-1 0;0 -1;1 2],[0;0;4]);
%    P_ = normalizedRep(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      13-September-2021 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    scal = diag(1./sqrt(sum(P.P.A.^2,2)));
    P = mptPolytope(scal * P.P.A, scal * P.P.b);

%------------- END OF CODE --------------