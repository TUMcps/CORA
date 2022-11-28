function S = simplex(Z)
% simplex - enclose a zonotope by a simplex
%
% Syntax:  
%    S = simplex(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    S - simplex represented as a mptPolytope object
%
% Example: 
%    Z = zonotope.generateRandom('Dimension', 2);
%    S = simplex(Z)
%    
%    figure; hold on; box on;
%    plot(Z);
%    plot(S, [1,2], 'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/interval, zonotope/mptPolytope

% Author:        Niklas Kochdumper
% Written:       31-May-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % construct an n-dimensional standard simplex with origin 0
    n = dim(Z);
    V = eye(n+1);
    B = gramSchmidt(ones(n+1,1));
    
    S = mptPolytope((B(:,2:end)'*V)');
    
    % scale the simplex so that it tightly encloses the zonotope    
    A = S.P.A;
    b = supremum(interval(A*Z));
    
    S = mptPolytope(A, b);

end

%------------- END OF CODE --------------