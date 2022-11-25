function res = conIntersect(Z,Y,R)
% conIntersect - Adds the constraint that the linear transformation R*Z of
%                a constrained zonotope Z has to intersect the constrained 
%                zonotope Y (see Eq. (13) in [1])
%
% Syntax:  
%    res = conIntersect(Z,Y,R)
%
% Inputs:
%    Z,Y - constrained zonotope object
%    R - transformation matrix for the linear transformation R*Z
%
% Outputs:
%    res - resulting constrained zonotope object
%
% Example: 
%    Z = conZonotope([0 1.5 -1.5 0.5;0 1 0.5 -1],[1 1 1],-1);
%    Y = conZonotope([6 2 0;-2 0 1]);
%    R = [2 1;-1 0];
%
%    res = conIntersect(Z,Y,R);
%
%    figure; hold on
%    plot(Z);
%    plot(res,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: and, linearSysDT/observe
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:        Niklas Kochdumper
% Written:       04-March-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------
       
    % get object properties
    cz = Z.Z(:,1); cy = Y.Z(:,1); Gz = Z.Z(:,2:end); Gy = Y.Z(:,2:end);
    
    % construct resulting set according to Eq. (13) in [1]
    c = cz;
    G = [Gz, zeros(size(Gz,1),size(Gy,2))];
    
    if ~isempty(Z.A)
       if ~isempty(Y.A)
           A = [blkdiag(Z.A,Y.A);[R*Gz -Gy]];
           b = [Z.b; Y.b; cy - R*cz];
       else
           A = [[Z.A,zeros(size(Z.A,1),size(Gy,2))];[R*Gz -Gy]];
           b = [Z.b; cy - R*cz];
       end
    else
       if ~isempty(Y.A)
           A = [[zeros(size(Y.A,1),size(Gz,2)) Y.A];[R*Gz -Gy]];
           b = [Y.b; cy - R*cz];
       else
           A = [R*Gz -Gy];
           b = cy - R*cz;
       end
    end

    res = conZonotope(c,G,A,b);
end

%------------- END OF CODE --------------