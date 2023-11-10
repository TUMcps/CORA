function [P,varargout] = polytope(Z,varargin)
% polytope - converts a zonotope object to a polytope object
%
% Syntax:
%    P = polytope(Z)
%    P = polytope(Z,'outer')
%
% Inputs:
%    Z - zonotope object
%    type - approximation:
%               'exact': based on Theorem 7 of [1]
%               'outer': ...
%    method - outer-approximation method for type = 'outer'
%               'tight': uses interval outer-approximation
%               'volume': uses volume
%
% Outputs:
%    P - polytope object
%    comb - (only method = 'exact') generator combinations corresponding to
%               the halfspaces
%    isDeg - (only method = 'exact') true/false whether polytope is
%               full-dimensional
%
% Example: 
%    Z = zonotope([1;-1],[3 2 -1; 1 -2 2]);
%    P = polytope(Z);
%
%    figure; hold on;
%    plot(P,[1,2],'r');
%    plot(Z,[1,2],'b');
%
% References:
%   [1] Althoff, M.; Stursberg, O. & Buss, M. Computing Reachable Sets
%       of Hybrid Systems Using a Combination of Zonotopes and Polytopes
%       Nonlinear Analysis: Hybrid Systems, 2010, 4, 233-249
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Matthias Althoff
% Written:       06-August-2018
% Last update:   10-November-2022 (MW, unify with other polytope functions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
[type,method] = setDefaultValues({'exact','tight'},varargin);

% check input arguments
inputArgsCheck({{Z,'att','zonotope'};
                {type,'str',{'exact','outer'}}; ...
                {method,'str',{'tight','volume'}}});


if strcmp(type,'exact')
    % note: this method was previously called 'polytope'

    % obtain number of generators, dimensions
    %Z = deleteAligned(Z);
    Z = compact_(Z,'zeros',eps);
    c = Z.c;
    G = Z.G;
    [n,nrGen] = size(G);
    
    isDeg = false;
    % other return values
    varargout{1} = [];
    
    if rank(G) >= n
        
        if n > 1
            % get number of possible facets
            comb = combinator(nrGen,n-1,'c');
            % bypass bug in combinator (rows with all-zeros?!)
            comb = comb(any(comb,2),:);
    
            % build C matrices for inequality constraint C*x < d
            C = zeros(length(comb(:,1)),n);
            for i=1:length(comb(:,1))
                indices=comb(i,:);
                Q=G(:,indices);
                v=ndimCross(Q);
                C(i,:)=v'/norm(v);
            end
    
            % remove NaN rows due to rank deficiency
            index = find(sum(isnan(C),2));
            if ~isempty(index)
                C(index,:) = [];
            end
        else
            C = 1;
            comb = [];
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
            
            % singular value decomposition
            [S,V,~] = svd(G);
            
            % state space transformation
            Z_ = S'*[c,G];
            
            % remove dimensions with all zeros
            ind = find(diag(V) <= 1e-12);
            ind_ = setdiff(1:size(V,1),ind);
            
            if ~isempty(ind)
                
                isDeg = true;
                
                % compute polytope in transformed space
                P = polytope(zonotope(Z_(ind_,:)));
                
                % transform back to original space
                A = [P.A,zeros(size(P.A,1),length(ind))]*S';
                b = P.b;
                
                % add equality constraint restricting polytope to null-space
                C = [A;S(:,ind)';-S(:,ind)'];
                d = [b;S(:,ind)'*c;-S(:,ind)'*c]; 
            end        
        end

        varargout{1} = comb;
        
    elseif nrGen == 0
        
        % generate equality constraint for the center vector
        E = eye(n);
        d_ = E * c;
        C = [E;-E];
        d = [d_;-d_];    
        
    else
        
        isDeg = true;
        
        % singular value decomposition
        [S,V,~] = svd(G);
        V = [V,zeros(n,n-nrGen)];
    
        % state space transformation
        Z_ = S'*[c,G];
    
        % remove dimensions with all zeros
        ind = find(diag(V) <= 1e-12);
        ind_ = setdiff(1:size(V,1),ind);
    
        if ~isempty(ind)
            
            % compute polytope in transformed space
            P = polytope(zonotope(Z_(ind_,:)));
    
            % transform back to original space
            A = [P.A,zeros(size(P.A,1),length(ind))]*S';
            b = P.b;
    
            % add equality constraint restricting polytope to null-space
            C = [A;S(:,ind)';-S(:,ind)'];
            d = [b;S(:,ind)'*c;-S(:,ind)'*c];
    
        end
    end
    
    % instantiate mpt Polytope
    P = polytope(C,d);
    % return values
    varargout{2} = isDeg;


elseif strcmp(type,'outer')
    % note: this method was previously called 'enclosingPolytope'
    
    if strcmp(method,'tight')
        % solution1 (axis-aligned):
        Zred = zonotope(interval(Z));
        P = polytope(Zred); 
        % solution 2 (method C):
%         Zred = reduce(Z,'methC',1,filterLength);
        Zred = reduce(Z,'pca');
        Zred = aux_repair(Zred,Z);
        Padd = polytope(Zred);
        % intersect results
        P = P & Padd;
    elseif strcmp(method,'volume')
        % solution 1 (method C):
%         Zred1 = reduce(Z,'methC',1,filterLength);
        Zred1 = reduce(Z,'pca');
        Zred1 = aux_repair(Zred1,Z);
        vol1 = volume(Zred1);
        % solution2 (axis-aligned):
        Zred2 = zonotope(interval(Z));
        Zred2 = aux_repair(Zred2,Z);
        vol2 = volume(Zred2);
    
        if vol1 < vol2
            P = polytope(Zred1);
        else 
            P = polytope(Zred2);
        end
    end

end

% remove redundant constraints too slow...
% P = compact_(P,'all',1e-9);

% Set polytope properties
P.bounded.val = true;
end


% Auxiliary functions -----------------------------------------------------

function Zrep = aux_repair(Z,Zorig)
% repair zonotope if there is no length in one dimenion

% get length of each dimension
len = 2*rad(interval(Z));

% find zero lengths
index = find(len==0);

if ~isempty(index)
    %construct zonotope to be added
    origLen = 2*rad(interval(Zorig));
    origCenter = center(interval(Zorig));
    
    %get Zmatrix
    Zmat = [Z.c,Z.G];
    
    for i=1:length(index)
        
        ind = index(i);
        %replace center
        Zmat(ind,1) = origCenter(ind);
        %replace generator value
        Zmat(ind,ind+1) = 0.5*origLen(ind);
    end

    % instantiate zonotopes
    Zrep = zonotope(Zmat);
else
    Zrep = Z;
end

end
    
% ------------------------------ END OF CODE ------------------------------
