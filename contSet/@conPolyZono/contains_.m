function [res,cert,scaling] = contains_(cPZ,S,method,tol,maxEval,certToggle,scalingToggle,varargin)
% contains_ - determines if a constrained polynomial zonotope contains a set
%    or a point
%
% Syntax:
%    [res,cert,scaling] = contains_(cPZ,S,method,tol,maxEval,certToggle,scalingToggle)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object or single point
%    method - method used for the containment check.
%       Currently, the only available options are 'exact' and 'approx'.
%    tol - tolerance for the containment check; the higher the
%       tolerance, the more likely it is that points near the boundary of
%       cPZ will be detected as lying in cPZ, which can be useful to
%       counteract errors originating from floating point errors.
%    maxEval - Currently has no effect
%    certToggle - if set to 'true', cert will be computed (see below),
%       otherwise cert will be set to NaN.
%    scalingToggle - if set to 'true', scaling will be computed (see
%       below), otherwise scaling will be set to inf.
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, S is
%           guaranteed to not be contained in cPZ, whereas if res=false and
%           cert=false, nothing can be deduced (S could still be
%           contained in cPZ).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(cPZ - cPZ.center) + cPZ.center contains S.
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    E = [1 0 3; 0 1 1; 0 0 0];
%    GI = [0.5;0]; 
%    A = [1 1 -1.5];
%    b = 0.5;
%    EC = [1 0 0;0 1 0;0 0 1];
%    cPZ = conPolyZono(c,G,E,A,b,EC);
% 
%    c = [0;0];
%    G = [1 -2 1; 2 3 1];
%    E = [1 0 2;0 1 1; 0 0 0];
%    A = [1 -1 0.5];
%    b = -0.5;
%    EC = [2 0 0;0 1 0; 0 0 1];
%    cPZ1 = conPolyZono(c,0.3*G,E,A,b,EC);
%    cPZ2 = conPolyZono(c,0.5*G,E,A,b,EC);
% 
%    p1 = [1;1];
%    p2 = [-1;3];
% 
%    contains(cPZ,p1,'approx')
%    contains(cPZ,p2,'approx')
%    contains(cPZ,cPZ1,'approx')
%    contains(cPZ,cPZ2,'approx')
% 
%    figure; hold on;
%    plot(cPZ,[1,2],'b','Splits',12);
%    plot(p1(1),p1(2),'.g','MarkerSize',20);
%    plot(p2(1),p2(2),'.r','MarkerSize',20);
%    plot(cPZ1,[1,2],'g','Splits',12);
%    plot(cPZ2,[1,2],'r','Splits',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains, interval/contains_, conZonotope/contains_

% Authors:       Niklas Kochdumper, Mark Wetzlinger, Adrian Kulmburg
% Written:       05-February-2020 
% Last update:   25-November-2022 (MW, rename 'contains')
%                16-January-2025 (AK, added scaling and cert)
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

if representsa(S, 'emptySet', tol)
    % Empty set is always contained
    res = true;
    cert = true;
    scaling = 0;
    return
elseif representsa(S, 'fullspace', tol)
    % Fullspace is never contained, since a cPZ is compact
    res = false;
    cert = true;
    scaling = inf;
    return
end

% The code is not yet ready to deal with scaling or cert
cert = false;
scaling = Inf;
if scalingToggle
    throw(CORAerror('CORA:notSupported',...
        "The computation of the scaling factor or cert " + ...
        "for constrained polynomial zonotopes is not yet implemented."));
end


% check user inputs 
if strcmp(method,'exact')
    throw(CORAerror('CORA:noExactAlg',cPZ,S));
end

% transform to equivalent higher-dimensional polynomial zonotope
m = size(cPZ.A,1);
c = [cPZ.c; -cPZ.b];
G = blkdiag(cPZ.G,cPZ.A);
E = [cPZ.E,cPZ.EC];

GI = cPZ.GI;
if ~isempty(GI)
    GI = [GI; zeros(m,size(GI,2))];
end

pZ = polyZonotope(c,G,GI,E,cPZ.id);

% increase dimension of the second set to match dim. of first set
if ~isempty(cPZ.A)
    if isnumeric(S)
        S = [S;ones(m,1)];
    else
        temp = sqrt(max(sum(cPZ.A.^2,1)))/100 * ones(m,1);
        S = cartProd_(S,interval(-temp,temp),'exact');
    end
end

% check containment using the function for higher-dimensional zonotopes
res = contains_(pZ,S,'approx',tol,maxEval,certToggle,scalingToggle);

% ------------------------------ END OF CODE ------------------------------
