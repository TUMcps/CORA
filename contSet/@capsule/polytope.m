function P = polytope(C,varargin)
% polytope - converts a capsule to a polytope object
%
% Syntax:
%    P = polytope(C)
%    P = polytope(C,method)
%
% Inputs:
%    C - capsule object
%    method - type of conversion: 'exact' (default), 'outer', 'inner'
%
% Outputs:
%    P - polytope object
%
% Example: 
%    C = capsule([1; 1], [0; 1], 0.5);
%    P = polytope(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       25-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default method: 'exact'
method = setDefaultValues({'exact'},varargin);

% check input arguments
inputArgsCheck({{C,'att','capsule'},...
                {method,'str',{'exact','outer','inner'}}});

% dimension
n = dim(C);

% empty case
if representsa_(C,'emptySet',0)
    P = polytope.empty(n);
    return
end

% exact conversion
if n == 1
    % 1D sets are intervals
    P = aux_1D(C);
    return
elseif all(withinTol(C.g,0,eps)) && withinTol(C.r,0,eps)
    % just a point
    P = polytope(C.c);
    return
elseif withinTol(C.r,0,eps)
    % lines in nD can also be exactly represented by polytopes
    P = aux_exactLineConversion(C,n);
    return
end

switch method
    case 'exact'
        % exact conversion only possible if the radius is 0 (see above)
        throw(CORAerror('CORA:noExactAlg'));

    case 'inner'
        throw(CORAerror('CORA:notSupported'));

    case 'outer'
        P = aux_supportFunction(C,n);

end

end


% Auxiliary functions -----------------------------------------------------

function P = aux_1D(C)

    % simple method for 1D case
    Vmin = C.c - abs(C.g) - C.r;
    Vmax = C.c + abs(C.g) + C.r;
    P = polytope([1;-1],[Vmax; -Vmin]);
    
    % set properties
    P.bounded.val = true;
    P.emptySet.val = false;
    P.minHRep.val = true;
    P.fullDim.val = ~withinTol(Vmin,Vmax,eps);
    if P.fullDim.val
        P.V.val = [Vmin, Vmax];
    else
        P.V.val = Vmin;
    end
    P.minVRep.val = true;

end

function P = aux_exactLineConversion(C,n)
% converts a line to a polytope

    % compute orthonormal basis
    B = gramSchmidt(C.g);

    % init polytope using two inequality constraints for the capsule
    % generator and equality generators for all orthogonal directions
    P = polytope([B(:,1)'; -B(:,1)'],[vecnorm(C.g);vecnorm(C.g)],...
                  B(:,2:end)',zeros(n-1,1));

    % shift by center
    P = P + C.c;

end

function P = aux_supportFunction(C,n)
% polytope enclosure via support function evaluations

    % compute orthonormal basis
    B = gramSchmidt(C.g);

    % enclose {C.g * beta | -1 <= beta <= 1} + {x | ||x||_2 <= r} by a
    % polytope via support function evaluations
    Z = zonotope(zeros(n,1),C.g);

    sF = zeros(2*n,1);
    % loop over all basis vectors
    for i=1:n
        % compute support function in positive and negative directions
        % (note that the support function of a ball is just the radius)
        sF(i) = supportFunc_(Z,B(:,i),'upper') + C.r;
        sF(n+i) = supportFunc_(Z,-B(:,i),'upper') + C.r;
    end

    % init polytope
    P = polytope([B'; -B'],sF);

    % shift polytope by center of capsule
    P = P + C.c;

end

% ------------------------------ END OF CODE ------------------------------
