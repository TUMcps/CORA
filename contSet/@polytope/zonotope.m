function Z = zonotope(P,varargin)
% zonotope - converts a polytope to a zonotope
%
% Syntax:
%    Z = zonotope(P)
%    Z = zonotope(P,mode)
%
% Inputs:
%    P - polytope object
%    mode - 'outer': outer approximation
%           'exact': exact conversion (not supported)
%           'inner': inner approximation
%
% Outputs:
%    Z - zonotope object
%
% Example: 
%    A = [1 0; -1 1; -1 -1]; b = [1;1;1];
%    P = polytope(A,b);
%    Z = zonotope(P);
%
%    figure; hold on
%    plot(P);
%    plot(Z,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       17-March-2023
% Last update:   08-January-2024 (MW, empty case)
%                13-March-2024 (TL, bug fix v rep given)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

mode = setDefaultValues({'outer'},varargin);
inputArgsCheck({{P,'att','polytope'},...
                {mode,'str',{'exact','outer','inner'}}});

% read out dimension
n = dim(P);

% empty set and fullspace
if representsa_(P,'emptySet',eps)
    Z = zonotope.empty(n);
    return
elseif representsa_(P,'fullspace',0) && n > 0
    % conversion of fullspace object not possible
    throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
        'can therefore not be converted into a zonotope.']));
end

% conversion requires constraints
constraints(P);

% different modes
switch mode
    case 'exact'
        throw(CORAerror('CORA:notSupported'));

    case 'inner'
        throw(CORAerror('CORA:notSupported'));

    case 'outer'
        Z = aux_zonotope_outer(P,n);
end

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_zonotope_outer(P,n)
% outer approximation of a polytope by a zonotope (parallelotope)

% compute maximum-volume ellipsoid to get "shape" of P and center
E = ellipsoid(P,'inner');
c = center(E);
Q = E.Q;
Q_r = sqrtm(Q);

% shift and transform polytope so that it is roughly a hyperbox
Mt = polytope(P);
Mt.A_.val = P.A_.val*Q_r;
Mt.b_.val = P.b_.val - P.A_.val*c;
if ~isempty(P.Ae_.val)
    Mt.Ae_.val = P.Ae_.val*Q_r;
    Mt.be_.val = P.be_.val - P.Ae_.val*c;
end
% also transform V rep if present
if P.isVRep.val
    Mt.V_.val = pinv(Q_r) * (P.V - c);
end

% compute enclosing hyperbox
sF_upper = zeros(n,1);
sF_lower = zeros(n,1);
I = eye(n);
for i=1:n
    sF_upper(i) = supportFunc(Mt,I(:,i),'upper');
    sF_lower(i) = supportFunc(Mt,I(:,i),'lower');
end
ct = 1/2*(sF_upper+sF_lower);
rt = 1/2*(sF_upper-sF_lower);

% init resulting zonotope
Z = Q_r*zonotope(ct,diag(rt)) + c;

end

% ------------------------------ END OF CODE ------------------------------
