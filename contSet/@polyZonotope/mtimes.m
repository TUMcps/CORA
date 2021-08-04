function pZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%          interval matrix with a polyZonotope
%
% Syntax:  
%    pZ = mtimes(matrix,pZ)
%
% Inputs:
%    matrix - numerical or interval matrix
%    pZ - polyZonotope object 
%
% Outputs:
%    pZ - polyZonotpe after multiplication of a matrix with a polyZonotope
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%    matrix = [1 2;-1 3];
%    intMatrix = interval([0.9 1.9;-1.1 2.9],[1.1 2.1;-0.9 3.1]);
%       
%    pZres = matrix*pZ;
%    pZresInt = intMatrix*pZ;
%   
%    figure
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    figure
%    plot(pZres,[1,2],'b','Filled',true,'EdgeColor','none');
%    figure
%    plot(pZresInt,[1,2],'g','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus, zonotope/mtimes

% Author:       Niklas Kochdumper
% Written:      25-June-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% find a poly zonotope object
if isa(factor1,'polyZonotope')
    
    pZ=factor1;
    matrix=factor2;

elseif isa(factor2,'polyZonotope')
    
    pZ=factor2;
    matrix=factor1;  
end

% numeric matrix
if isnumeric(matrix)
    
    pZ.c = matrix*pZ.c;
    
    if ~isempty(pZ.G)
        pZ.G = matrix*pZ.G;
    end
    
    if ~isempty(pZ.Grest)
        pZ.Grest = matrix*pZ.Grest;
    end

 
% interval matrix
elseif isa(matrix,'interval')
    
    m = center(matrix);
    r = rad(matrix);
    
    % calculate interval over-approximation
    inter = interval(pZ);
    s = abs(center(inter)) + rad(inter);
    
    % compute new polyZonotope
    pZ.c = m*pZ.c;
    if ~isempty(pZ.G)
        pZ.G = m*pZ.G;
    end
    if ~isempty(pZ.Grest)
        pZ.Grest = [m*pZ.Grest, diag(r*s)];
    else
        pZ.Grest = diag(r*s);
    end
    
    
% interval matrix 
elseif isa(matrix,'intervalMatrix')
    
    % get minimum and maximum
    M_min=infimum(matrix.int);
    M_max=supremum(matrix.int); 
    
    % get center and radius of interval matrix
    m=0.5*(M_max+M_min);
    r=0.5*(M_max-M_min);
    
    % calculate interval over-approximation
    inter = interval(pZ);
    s = abs(center(inter)) + rad(inter);
    
    % compute new polyZonotope
    pZ.c = m*pZ.c;
    if ~isempty(pZ.G)
        pZ.G = m*pZ.G;
    end
    if ~isempty(pZ.Grest)
        pZ.Grest = [m*pZ.Grest, diag(r*s)];
    else
        pZ.Grest = diag(r*s);
    end

    
% matrix zonotope 
else
    
    error('quadZonotope/mtimes: operation not yet implemented!')
end    




%-----------