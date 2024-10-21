function E_cell = priv_lplus(E_cell,L,mode)
% priv_lplus - Computes the Minkowski sum of a list of ellipsoids such that
%    the resulting over-approximation is tight in given directions
%
% Syntax:
%    E_cell = priv_lplus(E_cell,L,mode)
%
% Inputs:
%    E_cell - ellipsoid object (cell-array)
%    L  - unit directions
%    mode - type of approximation: 'inner', 'outer'
%
% Outputs:
%    E_cell - ellipsoid array after Minkowski addition
%             (its length is equal to the number of unit directions)
%
% Example: 
%    E1 = ellipsoid([3 -1; -1 1],[1;0]);
%    E2 = ellipsoid([5 1; 1 2],[1;-1]);
%    l = [1;0];
%    % E = priv_lplus({E1,E2},l);
%
% References:
%    [1] https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Victor Gassmann
% Written:       15-March-2019
% Last update:   05-July-2022 (VG, class array support)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% single ellipsoid
if length(E_cell) == 1
    return;
end

% reverse order to avoid changing array size each iteration
for i=size(L,2):-1:1
    E_cell{i} = aux_lplus_single(E_cell,L(:,i),mode);
end

end


% Auxiliary functions -----------------------------------------------------

function E_l = aux_lplus_single(E_cell,l,mode)
% see [1]
n = dim(E_cell{1});

switch mode 
    case 'outer'
        % outer-approximation
        q = zeros(n,1);
        c = 0;
        Q_ = zeros(n);
        for i=1:length(E_cell)
            q = q + E_cell{i}.q;
            si = sqrt(l'*E_cell{i}.Q*l);
            c = c + si;
            if ~withinTol(si,0,E_cell{i}.TOL)
                Q_ = Q_ + E_cell{i}.Q/si;
            end
        end
        Q = c*Q_;
        E_l = ellipsoid(Q,q);
    case 'inner'
        % inner-approximation
        x = sqrtm(E_cell{1}.Q)*l;
        q = zeros(n,1);
        Q = zeros(n);
        for i=1:length(E_cell)
            q = q + E_cell{i}.q;
            Qs = sqrtm(E_cell{i}.Q);
            Q = Q + vecalign(x,Qs*l)*Qs;
        end
        E_l = ellipsoid(Q'*Q,q);
end

end

% ------------------------------ END OF CODE ------------------------------
