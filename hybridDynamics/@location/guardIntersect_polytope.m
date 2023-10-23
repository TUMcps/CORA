function R = guardIntersect_polytope(loc,R,guard,options)
% guardIntersect_polytope - enclosure of guard intersections based on
%    using a combination of zonotopes and polytopes as described in [1]
%
% Syntax:
%    R = guardIntersect_polytope(loc,R,options)
%
% Inputs:
%    loc - location object
%    R - list of intersections between the reachable set and the guard
%    guard - guard set
%    options - struct containing the algorithm settings
%
% Outputs:
%    R - enclosure of the guard intersection
%
% References: 
%   [1] M. Althoff et al. "Computing Reachable Sets of Hybrid Systems Using 
%       a Combination of Zonotopes and Polytopes", 2009
%   [2] M. Althoff et al. "Zonotope bundles for the efficient computation 
%       of reachable sets", 2011

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       19-December-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % enclose all relevant reachable sets by polytopes
    for i = 1:length(R)
        R{i} = aux_conv2polytope(R{i});
    end

    % intersect the reachable sets with the guard set
    for i = 1:length(R)
        R{i} = and_(R{i},guard,'exact'); 
    end

    % compute vertices
    V = [];
    for i = 1:length(R)
        vert = vertices(R{i});
        V = [V,vert];
    end
    
    % enclose vertices with the methods described in Section V.A in [2]
    m = length(options.enclose);
    Z = cell(m,1);
    
    for i=1:m
        
        switch options.enclose{i}
            case 'box'
                % Box method as described in Section V.A.a) in [2]
                Z{i} = zonotope(interval.enclosePoints(V));
                
            case 'pca'
                % PCA mehtod as described in Section V.A.b) in [2]
                Z{i} = zonotope.enclosePoints(V,'stursberg');
                
            case 'flow'
                % Flow method as described in Section V.A.d) in [2]
                Z{i} = aux_flowEnclosure(loc.contDynamics,V,options); 
                
            otherwise
                throw(CORAerror('CORA:wrongFieldValue','options.enclose',...
                    {'box','pca','flow'}));
                
        end
    end
    
    % construct the enclosing zonotope or zonotope bundle object
    if length(Z) == 1
        R = Z{1}; 
    else
        R = zonoBundle(Z); 
    end 
end


% Auxiliary functions -----------------------------------------------------

function Z = aux_flowEnclosure(sys,V,options)
% Flow method for enclosing a set of vertices with a zonotope as described 
% in Section V.A.d) in [1]

    % compute flow at center of vertices
    c = mean(V,2);
    
    options.u = options.uTrans;
    fcnHan = getfcn(sys,options);
    dir = fcnHan(0,c);
    
    dir = dir./norm(dir);
    
    % get basis that is orthogonal to the flow direction
    B = gramSchmidt(dir);
    
    % compute box enclosure in transformed space 
    Z = B*zonotope(interval.enclosePoints(B'*V));

end

function P = aux_conv2polytope(R)
% compute an enclosing polytope for the given set

    % enclose set with zonotope
    if ~isa(R,'zonotope')
        R = zonotope(R); 
    end

    % reduce the given set with different methods
    Zred1 = reduce(R,'pca',1);
    Zred2 = zonotope(interval(R));
    
    % compute intersection of the two reduced sets
    P = and_(polytope(Zred1),polytope(Zred2),'exact');
end

% ------------------------------ END OF CODE ------------------------------
