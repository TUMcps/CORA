function P = reduce(P,method,order,varargin)
% reduce - reduces the order of a polynomial zonotope
%
% Syntax:
%    P = reduce(P,method,order,varargin)
%
% Inputs:
%    P - polytope object
%    option - str, reduction algorithm: 'rand'
%    order - order of reduced polytope
%
% Outputs:
%    pZ - reduced polytope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       08-December-2025 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(3,3)
admissibleMethods = {'rand'};
inputArgsCheck({ ...
    {P,'att','polytope'}, ...
    {method,'str',admissibleMethods}, ...
    {order,'att','numeric','scalar'}, ...
})

% compact polytope
P = compact(P);

% check if empty
if representsa(P,'emptySet')
    P = polytope.empty(dim(P)); return;
end

% reduce depending on method
switch method
    case 'rand'
        P = aux_rand(P,order);
        
    otherwise
        throw(CORAerror('CORA:wrongValue','second',admissibleMethods))
end

end


% Auxiliary functions -----------------------------------------------------

function P = aux_rand(P,order)
    % generate random directions and compute support functions

    n = dim(P);
    numDirs = 2*order*n;

    A = randn(numDirs,n);
    b = nan(numDirs,1);
    
    % compute support functions
    for i=1:size(A,1)
        b(i) = supportFunc_(P,A(i,:),'upper');
    end

    P = polytope(A,b);
    
end

% ------------------------------ END OF CODE ------------------------------
