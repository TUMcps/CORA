function [D_ksi,R] = priv_ksi_iterative(cZ,varargin)
% priv_ksi_iterative - determine the tightened domains for the zonotope
%    factors with an iterative method based on interval arithmetic
%
% Syntax:
%    Dksi = priv_ksi_iterative(cZ)
%    Dksi = priv_ksi_iterative(cZ,iter)
%
% Inputs:
%    cZ - conZonotope object
%    iter - number of performed iterations 
%
% Outputs:
%    D_ksi - new tightened domains for the zonotope factors (class: interval)
%    R - measure of how large the over-approximation would be if the
%        corresponding factor would be removed from the set representation
%        (needed for reduction of constrained zonotopes)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:       11-May-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
iter = setDefaultValues({1},varargin);

% check input arguments
inputArgsCheck({{cZ,'att','conZonotope'};
                {iter,'att','numeric','nonnan'}});

% initialization
[nrCon,nrGens] = size(cZ.A);
E = interval(-ones(nrGens,1), ones(nrGens,1));
R = interval(-Inf(nrGens,1), Inf(nrGens,1));

A = cZ.A;
b = cZ.b;
iA = A.^-1;

% loop over all iterations
for k = 1:iter
    % loop over all constraints
    for i = 1:nrCon
        % loop over all factors ksi
        for j = 1:nrGens
            if ~isinf(iA(i,j))
                % calculate new tightened domain for the current factor ksi
                temp = E;
                temp(j) = 0;
                dummy = iA(i,j) .* ( b(i) - A(i,:)*temp );
                
                % update domains
                R(j) = and_(R(j),dummy,'exact');
                E(j) = E(j) & R(j);
            end
        end
    end
end

D_ksi = E;

% ------------------------------ END OF CODE ------------------------------
