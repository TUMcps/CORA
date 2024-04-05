function cZ = minkDiff(cZ1,S,varargin)
% minkDiff - compute the Minkowski difference of two constrained zonotopes:
%         cZ1 - cZ2 = cZ <-> cZ + cZ2 \subseteq cZ1
%
% Syntax:
%    cZ = minkDiff(cZ1,S)
%    cZ = minkDiff(cZ1,S,method)
%
% Inputs:
%    cZ1 - conZonotope object
%    S - conZonotope object, contSet object, or numerical vector
%    method - type of computation
%           'exact' [1, Thm. 1]
%           'inner:vertices'
%           'inner:Vinod' [2, Alg. 1]
%
% Outputs:
%    cZ - conZonotope object after Minkowski difference
%
% Example: 
%    Z = [0 4 0 0;0 0 4 0];
%    A = [1 1 1.5]; b = -0.5;
%    cZ1 = conZonotope(Z,A,b);
%
%    Z = [1 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ2 = conZonotope(Z,A,b);
%
%    cZ = minkDiff(cZ1,cZ2);
%
%    figure; hold on;
%    plot(cZ1);
%    plot(cZ2,[1,2],'r');
%    plot(cZ,[1,2],'g');
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%    [2] A. Vinod, A. Weiss, S. Di Cairano, "Projection-free computation of
%        robust controllable sets with constrained zonotopes",
%        https://arxiv.org/pdf/2403.13730.pdf
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minkDiff

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       04-February-2021
% Last update:   09-November-2022 (MW, rename 'minkDiff')
%                02-April-2024 (MW, new method, refactor)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% different algorithms for different set representations
if isnumeric(S)
   cZ = cZ1 + (-S);
   return
end

% default method
method = setDefaultValues({'exact'},varargin);
admissibleMethods = {'exact','inner:vertices','inner:Vinod'};

switch method
    case 'exact'
        % compute exact result according to Lemma 1 in [1]
        V = vertices(S);
        cZ = cZ1 - V(:,1);
        for i = 2:size(V,2)
            cZ = and_(cZ,cZ1 - V(:,i),'exact'); 
        end

    case 'inner:vertices'
   
        % enclose subtrahend by zonotope
        Z = zonotope(S);
        
        % compute Minkowski difference according to Theorem 1 in [1]
        c = center(Z);
        G = generators(Z);
    
        cZ = cZ1 + (-c);
    
        for i = 1:size(G,2)
            cZ = and_(cZ + G(:,i),cZ + (-G(:,i)),'exact');
        end

    case 'inner:Vinod'
        % implements [2, Alg. 1]

        if ~isa(S,'zonotope') && ~isa(S,'interval')
            throw(CORAerror('CORA:notSupported'));
        end

        n = dim(cZ1);
        [M_C, N_C] = size(cZ1.A);
        c_S = center(S);
        S0 = S - c_S;
        Gamma = pinv([cZ1.G; cZ1.A]) * [eye(n); zeros(M_C,n)];

        S0_sF = zeros(N_C,1);
        for i=1:N_C
            basisvector = zeros(N_C,1);
            basisvector(i) = 1;
            S0_sF(i) = supportFunc_(S0,(basisvector'*Gamma)','upper');
            if S0_sF(i) > 1
                cZ = conZonotope.empty(n);
                return
            end
        end
        D = diag(ones(N_C,1) - S0_sF);
        cZ = conZonotope(cZ1.c - c_S, cZ1.G*D, cZ1.A*D, cZ1.b);

    otherwise
        throw(CORAerror('CORA:wrongValue','third',...
            "has to be " + strjoin(admissibleMethods,', ')));

end

% ------------------------------ END OF CODE ------------------------------
