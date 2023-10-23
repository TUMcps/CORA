function P_out = enclose(P1,varargin)
% enclose - computes a polytope that encloses a polytope and its linear 
%    transformation
%
% Syntax:
%    P_out = enclose(P1,P2)
%    P_out = enclose(P1,M,Pplus)
%
% Inputs:
%    P1 - first polytope object
%    P2 - second polytope object, satisfying P2 = (M * P1) + Pplus
%    M - matrix for the linear transformation
%    Pplus - polytope object added to the linear transformation
%
% Outputs:
%    P_out - polytope that encloses P1 and P2
%
% Example: 
%    A = [2 1; -2 2; -1 -2; 2 -1];
%    b = ones(4,1);
%    P1 = polytope(A,b);
%
%    M = [2 0; 0 1];
%    A = [1 0; 0 1; -1 -1];
%    b = [1; 1; -1];
%    Pplus = polytope(A,b);
%
%    P = enclose(P1,M,Pplus);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Authors:       Matthias Althoff, Viktor Kotsev
% Written:       02-February-2011
% Last update:   12-August-2016
%                31-May-2022
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin == 2
    P2 = varargin{1};
elseif nargin == 3
    P2 = (varargin{1}*P1) + varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% convex hull is equivalent to enclose (since result is convex)
P_out = convHull(P1,P2);

% ------------------------------ END OF CODE ------------------------------
