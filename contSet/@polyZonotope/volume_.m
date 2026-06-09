function vol = volume_(pZ,varargin)
% volume_ - Computes the volume of a polynomial zonotope
%
% Syntax:
%    vol = volume_(pZ)
%
% Inputs:
%    pZ - polyZonotope object
%
% Outputs:
%    vol - volume
%
% Example: 
%    pZ = polyZonotope([0;0],[2 1 2;0 2 2],[],[1 0 2;0 1 1]);
%    vol = volume(pZ)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/volume

% Authors:       Niklas Kochdumper
% Written:       07-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % redefine independent generators as new dependent generators
    if ~isempty(pZ.GI)
        pZ = polyZonotope(pZ.c,[pZ.G,pZ.GI],[], ...
                            blkdiag(pZ.E,eye(size(pZ.GI,2))));
    end

    % check if volume can be computed
    n = length(pZ.c); p = size(pZ.E,1);

    if p < n
        vol = 0; return;
    elseif p > n
        throw(CORAerror('CORA:notSupported',['Volume not supported ', ...
           'for polynomial zonotopes with more factors than dimensions']));
    end

    % compute Jacobian matrix of the polynomial zonotope function
    p = size(pZ.E,1);
    a = sym('a',[p,1]);

    f = fhandle(pZ);
    fsym = f(a);

    J = jacobian(fsym,a);

    % compute the determinant of the Jacobian matrix 
    d = det(J);

    % integrate over all variables
    for i = 1:p
        d = int(d,a(i));
        d = subs(d,a(i),1) - subs(d,a(i),-1);
    end

    % volume is exact if there are no overlaps (this is the case if
    % det(J(a)) ~= 0 forall a \in [-1,1]), and over-approximative otherwise
    vol = eval(d);

% ------------------------------ END OF CODE ------------------------------
