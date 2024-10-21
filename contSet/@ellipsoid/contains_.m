function res = contains_(E,S,mode,varargin)
% contains_ - determines if an ellipsoid contains a set or a point
%
% Syntax:
%    res = contains_(E,S)
%    res = contains_(E,S,mode)
%
% Inputs:
%    E - ellipsoid object 
%    S - contSet object or single point
%    mode - mode of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    E1 = ellipsoid([5 7;7 13],[1;2]);
%    E2 = ellipsoid(0.3*eye(2));
%    Z = zonotope([0 1 0;0 1 1]);
%
%    contains(E1,E2)
%    contains(E1,Z)
%
%    figure; hold on
%    plot(E1,[1,2],'b');
%    plot(E2,[1,2],'g');
%    plot(Z,[1,2],'r');
%
% References:
%    [1] Yildirim, E.A., 2006. On the minimum volume covering ellipsoid of
%        of ellipsoids. SIAM Journal on Optimization, 17(3), pp.621-641.     
%    [2] SDPT3: url: http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Victor Gassmann, Niklas Kochdumper, Adrian Kulmburg
% Written:       15-October-2019 
% Last update:   21-November-2019 (NK, extend to other sets)
%                09-March-2021 (included tolerance for q comparison)
%                17-March-2021 (error handling)
%                06-July-2021 (AK, merged containsPoint to in)
%                04-July-2022 (VG, handle class array cases)
%                25-November-2022 (MW, rename 'contains')
%                25-April-2023 (VG, add method for capsule)
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% containment check
if isnumeric(S)
    res = priv_containsPoint(E,S);
    return
end

if isa(S,'ellipsoid')
    res = priv_containsEllipsoid(E,S);
    return
end

if isa(S,'capsule')
    % check if balls at both ends of capsule are contained
    E1 = ellipsoid(S.r^2*eye(dim(S)), S.c+S.g);
    E2 = ellipsoid(S.r^2*eye(dim(S)), S.c-S.g);
    res = priv_containsEllipsoid(E,E1) && priv_containsEllipsoid(E,E2);
    return
end

if isa(S,'zonotope') && strcmp(mode,'approx')
    S = ellipsoid(S,'outer:norm_bnd');
    res = priv_containsEllipsoid(E,S);
    return
end

if strcmp(mode,'exact')
    if ismethod(S,'vertices_')
        % check if all vertices of the set are contained
        res = all(priv_containsPoint(E,vertices(S)));
        return
    end
    throw(CORAerror('CORA:noExactAlg',E,S));
end

% approx
if strcmp(mode,'approx')
    if ismethod(S,'zonotope')
        % zonotope over-approximation
        S = zonotope(S);
        res = contains_(E,S,mode,0);
        return
    end
end

throw(CORAerror('CORA:noops',E,S));

% ------------------------------ END OF CODE ------------------------------
