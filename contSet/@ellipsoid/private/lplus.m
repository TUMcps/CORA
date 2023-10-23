function E_L = lplus(E,L,mode)
% lplus - Computes the Minkowski sum of a list of ellipsoids such that the
%    resulting over-approximation is tight in given directions
%
% Syntax:
%    E_L = lplus(E,L,mode)
%
% Inputs:
%    E - ellipsoid object (array)
%    L  - unit directions
%    mode - type of approximation: 'inner', 'outer'
%
% Outputs:
%    E_L - Ellipsoid array after Minkowski addition (length(E_L)=size(L,2))
%
% Example: 
%    E1 = ellipsoid([3 -1; -1 1],[1;0]);
%    E2 = ellipsoid([5 1; 1 2],[1;-1]);
%    l = [1;0];
%    E = lplus([E1,E2],l);
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

if length(E)==1
    E_L = E;
    return;
end

% reverse order to avoid changing array size each iteration
for i=size(L,2):-1:1
    E_L(i) = aux_lplus_single(E,L(:,i),mode);
end

end


% Auxiliary functions -----------------------------------------------------

function E_l = aux_lplus_single(E,l,mode)
% see [1]
n = dim(E(1));
if strcmp(mode,'outer')
    % outer-approximation
    q = zeros(n,1);
    c = 0;
    Q_ = zeros(n);
    for i=1:length(E)
        q = q + E(i).q;
        si = sqrt(l'*E(i).Q*l);
        c = c + si;
        if ~withinTol(si,0,E(i).TOL)
            Q_ = Q_ + E(i).Q/si;
        end
    end
    Q = c*Q_;
    E_l = ellipsoid(Q,q);
else
    % inner-approximation
    x = sqrtm(E(1).Q)*l;
    q = zeros(n,1);
    Q = zeros(n);
    for i=1:length(E)
        q = q + E(i).q;
        Qs = sqrtm(E(i).Q);
        Q = Q + vecalign(x,Qs*l)*Qs;
    end
    E_l = ellipsoid(Q'*Q,q);
end

end

% ------------------------------ END OF CODE ------------------------------
