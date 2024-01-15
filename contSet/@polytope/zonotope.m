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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

mode = setDefaultValues({'outer'},varargin);
inputArgsCheck({{P,'att','polytope'},...
                {mode,'str',{'exact','outer','inner'}}});

% read out dimension
n = dim(P);

if representsa_(P,'emptySet',eps)
    Z = zonotope(zeros(n,0));
    return
end

if representsa_(P,'fullspace',0) && n > 0
    % conversion of fullspace object not possible
    throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
        'can therefore not be converted into a zonotope.']));
end

switch mode
    case 'exact'
        throw(CORAerror('CORA:notSupported'));

    case 'inner'
        throw(CORAerror('CORA:notSupported'));

    case 'outer'
        % compute maximum-volume ellipsoid to get "shape" of P and center
        E = ellipsoid(P,'inner');
        c = center(E);
        Q = E.Q;
        Q_r = sqrtm(Q);
        
        % shift and transform polytope so that it is roughly a hyperbox
        Mt = P;
        Mt.A = P.A*Q_r;
        Mt.b = P.b -P.A*c;
        if ~isempty(P.Ae)
            Mt.Ae = P.Ae*Q_r;
            Mt.be = P.be - M.Ae*c;
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
