function B = calcBasis(loc,R,guard,options,params)
% calcBasis - calculate orthogonal basis with the methods from [1]
%
% Syntax:
%    B = calcBasis(loc,R,guard,options,params)
%
% Inputs:
%    loc - object of class location
%    R - list of intersections between the reachable set and the guard
%    guard - guard set (class: constrained hyperplane)
%    options - struct containing the algorithm settings
%    params - model parameters for evaluating dynamic equation
%             (only options.enclose = 'flow')
%
% Outputs:
%    B - cell-array containing the calculated basis
%
% References: 
%   [1] M. Althoff et al. "Zonotope bundles for the efficient computation 
%       of reachable sets", 2011
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: location/guardIntersect

% Authors:       Niklas Kochdumper
% Written:       05-November-2018             
% Last update:   20-November-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % initialization
    sys = loc.contDynamics;
    m = length(options.enclose);
    B = cell(m,1);
    
    % loop over all selected methods
    for i = 1:m
        
        switch options.enclose{i}
           
            % Box method as described in Section V.A.a) in [1]
            case 'box'
                B{i} = eye(sys.nrOfDims);
                
                
            % PCA mehtod as described in Section V.A.b) in [1]
            case 'pca'
                
                % concatenate all generators
                G = aux_extractGenerators(R);
                
                % project the generators onto the hyperplane
                if isa(guard,'polytope') && representsa_(guard,'conHyperplane',1e-12)
                    Z = zonotope([zeros(sys.nrOfDims,1),G]); 
                    Z_ = projectOnHyperplane(Z,guard);
                    G = generators(Z_);
                end
                
                % limit maximum number of generators to 1000
                if size(G,2) > 500
                    [~,ind] = sort(sum(G.^2,1),'descend');
                    G = G(:,ind(1:500));
                end
                
                % calcualte an orthogonal transformation matrix using PCA
                [B{i},~,~] = svd([-G,G]); 
                
                
            % Flow method as described in Section V.A.d) in [1]
            case 'flow'
                
                % compute projected center of the union of all sets
                c = aux_extractCenter(R);
                
                if isa(guard,'polytope') && representsa_(guard,'conHyperplane',1e-12)
                   Z = zonotope(c); 
                   Z_ = projectOnHyperplane(Z,guard);
                   c = center(Z_);
                end
                
                % compute flow at the center
                params.u = params.uTrans;
                fcnHan = getfcn(sys,params);
                dir = fcnHan(0,c);

                dir = dir./norm(dir);

                % get basis that is orthogonal to the flow direction
                B{i} = gramSchmidt(dir);
            
            otherwise
                throw(CORAerror('CORA:wrongFieldValue','options.enclose',...
                    {'box','pca','flow'}));
        end        
    end
end


% Auxiliary functions -----------------------------------------------------

% TODO: move to global!

function G = aux_extractGenerators(R)
% extract all generator vectors from the sets that are over-approximated

    G = [];

    % loop over all sets
    for j = 1:length(R)
        
        % different types of sets
        if isa(R{j},'zonoBundle')
            for i = 1:R{j}.parallelSets
                G = [G,R{j}.Z{i}.G];
            end
        elseif isa(R{j},'polyZonotope')
            temp = zonotope(R{j});
            G = [G,temp.G];
        elseif isa(R{j},'zonotope') || isa(R{j},'conZonotope')
            G = [G,R{j}.G];
        else % ?
            G = [G,R{j}.Z(:,2:end)];
        end
    end
end

function c = aux_extractCenter(R)
% extract all generator vectors from the sets that are over-approximated

    c = [];

    % loop over all sets
    for j = 1:length(R)

        % different types of sets
        if isa(R{j},'zonoBundle')
            for i = 1:R{j}.parallelSets
                c = [c,R{j}.Z{i}.c];
            end
        else
            c = [c,center(R{j})];
        end
    end
    
    % compute center of all centers
    c = mean(c,2);
end

% ------------------------------ END OF CODE ------------------------------
