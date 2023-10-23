function E = plusEllipsoidOA_halder(E)
% plusEllipsoidOA_halder - Computes an approximation to the smallest-volume
%    outer-approximation of the Minkowski sum of ellipsoids
%
% Syntax:
%    E = plusEllipsoidOA_halder(E)
%
% Inputs:
%    E - array of ellipsoid objects
%
% Outputs:
%    E - ellipsoid object
%
% References:
%    [1] A. Halder. "On the parameterized computation of minimum volume
%        outer ellipsoid of Minkowski sum of ellipsoids.", CDC 2018
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Victor Gassmann
% Written:       25-May-2022
% Last update:   05-July-2022 (VG, class array support)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% no error checking, already done in parent function

% number of ellipsoids in cell array
N = length(E);

% read values from first ellipsoid
Q = E(1).Q;
q = E(1).q;
TOL = E(1).TOL;

% loop over ellipsoids
for i=2:N
    Q1 = Q;
    Q2 = E(i).Q;
    l_i = eig(Q1\Q2);
    
    beta = 0.5;
    while true
        beta_new = (sum(1./(1+beta*l_i))/sum(l_i./(1+beta*l_i)))^(1/2);
        if withinTol(beta,beta_new,TOL)
            break;
        end
        beta = beta_new;
    end

    if beta_new<0 || withinTol(beta_new,0,TOL)
        throw(CORAerror('CORA:specialError','Value almost zero!'));
    end

    Q = (1+1/beta_new)*Q1+(1+beta_new)*Q2;
    q = q + E(i).q;
end

% instantiate resulting ellipsoid
E = ellipsoid(Q,q);

% ------------------------------ END OF CODE ------------------------------
