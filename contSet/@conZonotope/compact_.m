function cZ = compact_(cZ,method,tol,varargin)
% compact_ - returns equal constrained zonotope in minimal representation
%    1) delete all constraints of the form 0 * beta = 0 as they are
%    trivially true for all values of beta; additionally, remove all
%    generators which have zero-length and no corresponding entries in the
%    constraint matrix
%    note that there are also other constraints which might be true for all
%    values of beta (-1 to 1), but this would require much more
%    computational effort and is thus omitted in this function
%
% Syntax:
%    cZ = compact_(cZ)
%
% Inputs:
%    cZ - conZonotope object
%    method - method for redundancy removal
%             'zeros' (default): delete all constraints of the form
%                  0 * beta = 0 as they are trivially true for all values
%                  of beta; additionally, remove all generators which have
%                  zero-length and no corresponding entries in the
%                  constraint matrix (note that there are also other
%                  constraints which might be true for all values of beta
%                  (-1 to 1), but this would require much more
%                  computational effort and is thus omitted here)
%    tol - tolerance
%
% Outputs:
%    cZ - conZonotope object
%
% Example: 
%    Z = [0 1 0 0 1;0 1 0 2 -1];
%    A = [-2 0 1 -1; 0 0 0 0]; b = [2;0];
%    cZ = conZonotope(Z,A,b);
%    
%    cZ_ = compact(cZ,'zeros');
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    plot(cZ_,[1,2],'b--');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/compact, zonotope/compact_
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       29-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

switch method
    case 'all'
        cZ = compact_(cZ,'zeros',tol,varargin{:});

    case 'zeros'
        
        % constraints that are trivially true: 0 * beta = 0
        idx = ~any([cZ.A,cZ.b],2);
        cZ.A(idx,:) = []; cZ.b(idx) = [];
        
        % zero-length generators (corresponding columns in constraint matrix A need
        % to be all-zero as well)
        idx = ~any([cZ.G;cZ.A],1);
        % remove zero-length generators and corresponding constraints
        cZ.G(:,idx) = [];
        if ~isempty(cZ.A)
            cZ.A(:,idx) = [];
        end

end

% ------------------------------ END OF CODE ------------------------------
