function cPZ = conPolyZono(E)
% conPolyZono - Convert an ellipsoid to a conPolyZono object
%
% Syntax:  
%    cPZ = conPolyZono(E)
%
% Inputs:
%    E - ellipsoid object
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    E = ellipsoid([4 2;2 4],[1;1]);
%    cPZ = conPolyZono(E);
% 
%    figure; hold on;
%    plot(E,[1,2],'r','Filled',true,'EdgeColor','none');
%    xlim([-2,4]); ylim([-2,4]);
% 
%    figure; hold on;
%    plot(cPZ,[1,2],'b','Splits',20,'Filled',true,'EdgeColor','none');
%    xlim([-2,4]); ylim([-2,4]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Author:       Niklas Kochdumper
% Written:      12-August-2019 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % eigenvalue decomposition of the ellipsoid matrix
    [V,D] = eig(E.Q);
    
    % construct constrained polynomial zonotopes
    n = E.dim;
    
    c = E.q;
    G = V*sqrt(D);
    expMat = [eye(n);zeros(1,n)];
    A = [-0.5 ones(1,n)];
    b = 0.5;
    expMat_ = [[zeros(n,1);1],[2*eye(n); zeros(1,n)]];
    id = (1:n+1)';
    
    cPZ = conPolyZono(c,G,expMat,A,b,expMat_,[],id);
    
end
    
%------------- END OF CODE --------------