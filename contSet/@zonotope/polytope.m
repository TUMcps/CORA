function [P,comb,isDeg] = polytope(Z, varargin)
% polytope - converts a zonotope from a G- to a H-representation
%    This function is implemented based on Theorem 7 of [1].
%
% Syntax:  
%    [P,comb,isDeg] = polytope(Z)
%    [P,comb,isDeg] = polytope(Z,type)
%
% Inputs:
%    Z - zonotope object
%    type - type of polytope returned
%               - 'mpt' (default)
%               - 'ppl'
%
% Outputs:
%    P - polytope object
%    comb - generator combinations corresponding to the halfspaces
%    isDeg - flag spacifying if polytope is full-dimensional (0) or not (1)
%
% Example: 
%    zono = zonotope.generateRandom(2,[],5);
%    poly = polytope(zono);
%
%    figure; hold on;
%    plot(poly,[1,2],'r');
%
%    figure; hold on;
%    plot(zono,[1,2],'b');
%
% References:
%   [1] Althoff, M.; Stursberg, O. & Buss, M. Computing Reachable Sets
%       of Hybrid Systems Using a Combination of Zonotopes and Polytopes
%       Nonlinear Analysis: Hybrid Systems, 2010, 4, 233-249
% 
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  vertices

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      30-September-2008
% Last update:  26-February-2009
%               05-July-2010
%               20-October-2010
%               03-December-2010
%               02-September-2011 (delete aligned added)
%               12-February-2012
%               13-October-2014
%               20-March-2015
%               02-June-2015
%               11-June-2015
%               31-August-2015
%               12-August-2016
%               30-September-2016 (one-dimensional case added)
%               22-January-2019 (NK, non full-dimensional case added)
%               28-October-2019 (NK, catch not full-dimensional case)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin == 1
    options.polytopeType = 'mpt';
elseif nargin == 2
    options = varargin{1};
end

% obtain number of generators, dimensions
%Z=deleteAligned(Z);
Z = deleteZeros(Z);
c = center(Z);
G = generators(Z);
[n,nrGen] = size(G);

isDeg = 0;

if nrGen >= n
    
    if n > 1
        % get number of possible facets
        comb = combinator(nrGen,n-1,'c');
        % bypass bug in combinator (rows with all-zeros?!)
        comb = comb(any(comb,2),:);

        % build C matrices for inequality constraint C*x < d
        C=[];
        for i=1:length(comb(:,1))
            indices=comb(i,:);
            Q=G(:,indices);
            v=ndimCross(Q);
            C(end+1,:)=v'/norm(v);
        end

        % remove NaN rows due to rank deficiency
        index = find(sum(isnan(C),2));
        if ~isempty(index)
            C(index,:) = [];
        end
    else
        C = 1;
    end
    
    % build d vector and determine delta d
    deltaD = zeros(length(C(:,1)),1);
    for iGen = 1:nrGen
        deltaD = deltaD+abs(C*G(:,iGen));
    end 

    %compute dPos, dNeg
    dPos = C*c + deltaD;
    dNeg = -C*c + deltaD;
     
    % construct the overall inequality constraints
    C = [C;-C];
    d = [dPos;dNeg];
    
    % catch the case where the zonotope is not full-dimensional
    temp = min([sum(abs(C - C(1,:)),2),sum(abs(C + C(1,:)),2)],[],2);
    
    if n > 1 && (isempty(C) || all(all(isnan(C))) || ...
                   all(temp < 1e-12) || any(max(abs(C),[],1) < 1e-12))
        
        % singluar value decomposition
        [S,V,~] = svd(G);
        
        % state space transformation
        Z_ = S'*[c,G];
        
        % remove dimensions with all zeros
        ind = find(diag(V) <= 1e-12);
        ind_ = setdiff(1:size(V,1),ind);
        
        if ~isempty(ind)
            
            isDeg = 1;
            
            % compute polytope in transformed space
            poly = polytope(zonotope(Z_(ind_,:)));
            P = get(poly,'P');
            
            % transform back to original space
            A = [P.A,zeros(size(P.A,1),length(ind))]*S';
            b = P.b;
            
            % add equality constraint restricting polytope to null-space
            C = [A;S(:,ind)';-S(:,ind)'];
            d = [b;S(:,ind)'*c;-S(:,ind)'*c]; 
        end        
    end        
    
elseif nrGen == 0
    
    % generate equality constraint for the center vector
    E = eye(n);
    d_ = E * c;
    C = [E;-E];
    d = [d_;-d_];    
    
else
    
    isDeg = 1;
    
    % singluar value decomposition
    [S,V,~] = svd(G);
    V = [V,zeros(n,n-nrGen)];

    % state space transformation
    Z_ = S'*[c,G];

    % remove dimensions with all zeros
    ind = find(diag(V) <= 1e-12);
    ind_ = setdiff(1:size(V,1),ind);

    if ~isempty(ind)
        
        % compute polytope in transformed space
        poly = polytope(zonotope(Z_(ind_,:)));
        P = get(poly,'P');

        % transform back to original space
        A = [P.A,zeros(size(P.A,1),length(ind))]*S';
        b = P.b;

        % add equality constraint restricting polytope to null-space
        C = [A;S(:,ind)';-S(:,ind)'];
        d = [b;S(:,ind)'*c;-S(:,ind)'*c];

    end
end

% convert to mpt or ppl Polytope
if isfield(options,'polytopeType') && strcmp(options.polytopeType,'ppl')
    P = pplPolytope(C,d);
else
    P = mptPolytope(C,d);
    %P = removeRedundancies(P,'aligned');
end


%------------- END OF CODE --------------
