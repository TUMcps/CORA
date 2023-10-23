function pZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%    interval matrix with a polynomial zonotope
%
% Syntax:
%    pZ = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or interval matrix
%    factor2 - polyZonotope object 
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
%    figure; hold on;
%    plot(pZ,[1,2],'r');
%    plot(pZres,[1,2],'b');
%    plot(pZresInt,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus, zonotope/mtimes

% Authors:       Niklas Kochdumper
% Written:       25-June-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% find the polyZonotope object
[pZ,matrix] = findClassArg(factor1,factor2,'polyZonotope');


try

    % numeric matrix
    if isnumeric(matrix)
        
        pZ.c = matrix*pZ.c;
        
        if ~isempty(pZ.G)
            pZ.G = matrix*pZ.G;
        end
        
        if ~isempty(pZ.GI)
            pZ.GI = matrix*pZ.GI;
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
        if ~isempty(pZ.GI)
            pZ.GI = [m*pZ.GI, diag(r*s)];
        else
            pZ.GI = diag(r*s);
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
        if ~isempty(pZ.GI)
            pZ.GI = [m*pZ.GI, diag(r*s)];
        else
            pZ.GI = diag(r*s);
        end
    
        
    % matrix zonotope 
    else
    
        throw(CORAerror('CORA:noops',factor1,factor2));
    
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check for empty sets
    if representsa_(pZ,'emptySet',eps)
        return
    end

    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);

    % other error...
    rethrow(ME);

end

% ------------------------------ END OF CODE ------------------------------
