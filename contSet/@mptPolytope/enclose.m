function P = enclose(P,varargin)
% enclose - encloses a polytope and its affine transformation
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * (M x1 + x2) | x1 \in P, x2 \in P2, a \in [0,1] }
%    where P2 = M*P + Pplus
%
% Syntax:  
%    P = enclose(P1,P2)
%    P = enclose(P1,M,Pplus)
%
% Inputs:
%    P1 - mptPolytope object
%    P2 - mptPolytope object
%    M - matrix for the linear transformation
%    Pplus - mptPolytope object added to the linear transformation
%
% Outputs:
%    obj - mptPolytope object that encloses P1 and P2
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Author:       Matthias Althoff
% Written:      02-February-2011
% Last update:  12-August-2016
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin == 2
    P2 = varargin{1};
elseif nargin == 3
    P2 = (varargin{1}*P) + varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% call MPT-toolbox method to over-approximate with convex hull
try %MPT3
    P.P = convexHull(P.P, P2.P);
catch %MPT2
    %compute convex hull
    P.P = hull([P.P P2.P]);
end

%------------- END OF CODE --------------