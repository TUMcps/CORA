function E_cell = lplus(E_c,L,mode)
% lplus - Computes the Minkowski sum of all ellipsoids contained in E_c such that resulting
% overapproximation given by E is tight in directions L
%
% Syntax:
%    E_cell = lplus(E_c,l)
%
% Inputs:
%    E_c - cell array containing ellipsoids
%    L  - unit directions
%    mode-"i": inner approx; "o"(default): outer approx
%
% Outputs:
%    E - Ellipsoid after minkowski addition
%
% Example: 
%    E1=ellipsoid.generateRandom(2,false);
%    E2=ellipsoid.generateRandom(2,false);
%    l =[1;0];
%    E =lplus({E1,E2},l);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus
%
% References:
%    [1] https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%
% Author:       Victor Gassmann
% Written:      15-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


if length(E_c)==1
    E_cell = E_c;
    return;
end

E_cell = cell(length(E_c),1);
for i=1:size(L,2)
    E_cell{i} = lplus_single(E_c,L(:,i),mode);
end

end


%--- helper
function E = lplus_single(E_c,l,mode)
n = length(E_c{1}.q);
if strcmp(mode,'o')
    % outer approximation
    q = zeros(n,1);
    c = 0;
    Q_ = zeros(n);
    for i=1:length(E_c)
        q = q + E_c{i}.q;
        si = sqrt(l'*E_c{i}.Q*l);
        c = c + si;
        if ~withinTol(si,0,E_c{i}.TOL)
            Q_ = Q_ + E_c{i}.Q/si;
        end
    end
    Q = c*Q_;
    E = ellipsoid(Q,q);
else
    % inner approximation
    x = sqrtm(E_c{1}.Q)*l;
    q = zeros(n,1);
    Q = zeros(n);
    for i=1:length(E_c)
        q = q + E_c{i}.q;
        Qs = sqrtm(E_c{i}.Q);
        Q = Q + vecalign(x,Qs*l)*Qs;
    end
    E = ellipsoid(Q'*Q,q);
end
end
%------------- END OF CODE --------------