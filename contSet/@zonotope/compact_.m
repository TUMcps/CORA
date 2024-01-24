function Z = compact_(Z,method,tol,varargin)
% compact_ - returns equal zonotope in minimal representation
%
% Syntax:
%    Z = compact_(Z)
%    Z = compact_(Z,method)
%    Z = compact_(Z,method,tol)
%
% Inputs:
%    Z - zonotope object
%    method - method for redundancy removal
%             'zeros' (default): removes all generators from a zonotope
%                  with zeros in all dimensions, so that every generator
%                  the resulting zonotope has at least one non-zero entry
%             'aligned': combines aligned generators to a single generator;
%                  a tolerance is used to determine alignment, so this
%                  function does not necessarily return an
%                  over-approximation of the original zonotope---for this,
%                  use zonotope/reduce instead
%             'all' - all methods in succession
%    tol - tolerance
%
% Outputs:
%    Z - zonotope object
%
% Example:
%    Z1 = zonotope([0;0],[1 0 -2 0 3 4; 0 0 1 0 -2 1]);
%    Z2 = compact(Z1);
%    
%    plot(Z1); hold on;
%    plot(Z2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/compact

% Authors:       Mark Wetzlinger, Matthias Althoff
% Written:       15-January-2009
% Last update:   27-August-2019
% Last revision: 29-July-2023 (MW, merged from deleteZeros/deleteAligned)

% ------------------------------ BEGIN CODE -------------------------------

switch method
    case 'zeros'
        Z = aux_deleteZeros(Z,tol);
    case 'aligned'
        Z = aux_deleteAligned(Z,tol);
    case 'all'
        Z = aux_deleteZeros(Z,tol);
        Z = aux_deleteAligned(Z,tol);
end

end


% Auxiliary functions -----------------------------------------------------

function Z = aux_deleteZeros(Z,tol)
    
    % filter zero generators
    Z.G = nonzeroFilter(Z.G,tol);

end

function Z = aux_deleteAligned(Z,tol)

    % extract generator matrix
    G = Z.G;
    
    % delete zero-generators
    G = nonzeroFilter(G);
    
    % normalize generators
    G_norm = G./vecnorm(G);
    
    % tolerance for alignment
    tol = 1-tol;
    
    % find equal generators
    i = 1;
    while i < length(G(1,:))
        G_act = G_norm(:,i);
        ind = find(abs(G_act'*G_norm(:,(i+1):end)) > tol);
        if ~isempty(ind)
            ind = ind+i;
            for iAdd = 1:length(ind)
                %add generators
                G(:,i) = G(:,i) + sign(G_act'*G_norm(:,ind(iAdd)))*G(:,ind(iAdd));
            end
            for iAdd = 1:length(ind)
                %remove generators
                G(:,ind(iAdd)) = [];
                G_norm(:,ind(iAdd)) = [];
                %change ind to correct for cancellation
                ind = ind - 1;
            end
        end  
        % increase i
        i = i + 1;
    end
    
    Z.G = G;

end

% ------------------------------ END OF CODE ------------------------------
