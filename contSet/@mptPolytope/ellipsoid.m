function E = ellipsoid(obj,mode)
% ellipsoid - Overapproximates a mptPolytope by an ellipsoid
%
% Syntax:  
%    E = ellipsoid(obj,comptype)
%
% Inputs:
%    Z       - zonotope object
%    mode    - (Optional) 'i' (inner approx); 'o' (outer approx)
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    Z = zonotope(rand(2,5));
%    E = ellipsoid(Z);%same as ellipsoid(Z,'o:norm:bnd')
%    plot(Z);
%    hold on
%    plot(E);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      15-March-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
if ~exist('mode','var')
    mode = 'o';
end

V = vertices(obj);
% remove bias
m = mean(V,2);
V = V-m;

if strcmp(mode,'o')
    E = ellipsoid.enclosePoints(V);
else
    % first, check whether obj is degenerate
    if isempty(V)
        E = ellipsoid;
        return;
    end
    n = size(V,1);
    [T,S,~] = svd(V);
    nt = rank(ellipsoid(S(:,1:n)));
    V = T'*V;
    n_d = n-nt;
    r_d = zeros(n_d,1);
    if nt<n
        % save rest
        r_d = 1/2*(max(V(nt+1:end,:),[],2)-min(V(nt+1:end,:),[],2));
        % remove zeros if degenerate
        V(nt+1:end,:) = [];
    end
    obj = mptPolytope(V');
    A = obj.P.A;
    b = obj.P.b;
    % normalize to prevent numerical issues
    fac = 1./sqrt(sum(A.^2,2));
    A = 1./sqrt(sum(A.^2,2)).*A;
    b = fac.*b;
    c = sdpvar(nt,1);
    B = sdpvar(nt);
    F = cone([b-A*c,A*B]');
    options = sdpsettings;
    options.verbose = 0;
    if exist('sdpt3','file')
        options.solver = 'sdpt3';
    end
    %solve optimization problem
    diagnostics = optimize(F, -logdet(B), options);
    if diagnostics.problem ~= 0
        throw(errOptNumIssue());
    end
    E_ext = ellipsoid(value(B)^2, value(c));
    for i=1:n_d
        E_ext = cartProd(E_ext,ellipsoid(r_d(i)^2));
    end
    E = T*E_ext;
end
% add back center
E = E + m;
%------------- END OF CODE --------------