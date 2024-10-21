function E = priv_plusEllipsoidOA_halder(E_cell)
% priv_plusEllipsoidOA_halder - computes an approximation to the
%    smallest-volume outer approximation of the Minkowski sum of ellipsoids
%
% Syntax:
%    E = priv_plusEllipsoidOA_halder(E_cell)
%
% Inputs:
%    E_cell - cell-array of ellipsoid objects
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
% See also: ellipsoid/plus, priv_plusEllipsoidOA

% Authors:       Victor Gassmann
% Written:       25-May-2022
% Last update:   05-July-2022 (VG, class array support)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% number of ellipsoids in cell array
N = length(E_cell);

% read values from first ellipsoid
Q = E_cell{1}.Q;
q = E_cell{1}.q;
TOL = E_cell{1}.TOL;

% loop over ellipsoids
for i=2:N
    Q1 = Q;
    Q2 = E_cell{i}.Q;
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
    q = q + E_cell{i}.q;
end

% instantiate resulting ellipsoid
E = ellipsoid(Q,q);

% ------------------------------ END OF CODE ------------------------------
