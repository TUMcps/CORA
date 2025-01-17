function cZ = conZonotope(P,varargin)
% conZonotope - converts a polytope object into a constrained zonotope;
%    implementation according to Theorem 1 from [1]
%
% Syntax:
%    cZ = conZonotope(P)
%    cZ = conZonotope(P,method)
%    cZ = conZonotope(P,method,box)
%
% Inputs:
%    P - polytope object
%    method - (optional) conversion method
%             - 'exact:supportFunc': using support functions (default)
%             - 'exact:vertices': using vertex enumeration
%    box - (optional) box outer approximation of P
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    A = [-1 0;0 -1;1 1];
%    b = [1;1;1];
%    P = polytope(A,b);
%    cZ = conZonotope(P);
%
%    figure; hold on
%    plot(cZ);
%    plot(P,[1,2],'r--');
%
% References: 
%    [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%        estimation and fault detection"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       13-May-2018
% Last update:   28-April-2019 (MA, code shortened)
%                13-December-2022 (MW, add support-function based method)
%                14-July-2023 (MW, add support for empty polytopes)
%                27-July-2023 (MW, incorporate equality constraints)
%                03-January-2024 (MW, speed up supportFunc method, handle unbounded)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,3);

% dimension
n = dim(P);

% set default method
[method,B] = setDefaultValues({'exact:supportFunc',interval.empty(n)},varargin);

% check input arguments
inputArgsCheck({{P,'att','polytope'}, ...
                {method,'str',{'exact:vertices','exact:supportFunc'}}, ...
                {B,'att','interval'}});

% fullspace cannot be represented as a conZonotope
if representsa_(P,'fullspace',0)
    % conversion of fullspace object not possible
    throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
        'can therefore not be converted into a constrained zonotope.']));
end

% vertex representation given and polytope is empty
if P.isVRep.val && isempty(P.V_.val)
    cZ = conZonotope.empty(n);
    return
end

% choose method
switch method
    case 'exact:vertices'
        % calculate the vertices of the polytope (also check for unboundedness)
        try
            V = vertices(P);
        catch ME
            if ~isempty(P.bounded.val) && ~P.bounded.val
                 throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
                    'can therefore not be converted into a constrained zonotope.']));
            end
            rethrow(ME);
        end
        
        % no vertices -> empty set
        if isempty(V)
            cZ = conZonotope.empty(n);
            return
        end

        % ensure that constraints are there
        [A,b,Ae,be] = constraints(P);
        
        % conversion
        [c_,G_,A_,b_] = priv_conZonotope_vertices(A,b,Ae,be,V);

    case 'exact:supportFunc'
        % compute bounding box (which also computes set properties)
        if nargin < 3
            B = interval(P);
        end

        % check if P is known to be unbounded or empty
        if ~isempty(P.bounded.val) && ~P.bounded.val
            throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
                'can therefore not be converted into a constrained zonotope.']));
        elseif ~isempty(P.emptySet.val) && P.emptySet.val
            cZ = conZonotope.empty(n);
            return
        end

        % ensure that constraints are there
        [A,b,Ae,be] = constraints(P);
        
        % conversion
        [c_,G_,A_,b_,empty] = priv_conZonotope_supportFunc(A,b,Ae,be,B);
        if empty
            cZ = conZonotope.empty(n);
            return
        end

end

% init constained zonotope
cZ = conZonotope(c_,G_,A_,b_);

% ------------------------------ END OF CODE ------------------------------
