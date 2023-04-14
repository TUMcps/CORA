function cPZ = mtimes(factor1,factor2)
% mtimes - Overloaded '*' operator for the multiplication of a matrix (or
%    interval matrix) with a constrained polynomial zonotope
%
% Syntax:  
%    cPZ = mtimes(factor1,factor2)
%
% Inputs:
%    factor1 - numerical matrix or interval matrix
%    factor2 - conPolyZono object 
%
% Outputs:
%    cPZ - conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [2 0 2; 0 2 2];
%    expMat = [1 0 1; 0 1 1; 0 0 0];
%    A = [2 2 4 -4];
%    b = 0;
%    expMat_ = [1 0 1 0; 0 1 1 0; 0 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
% 
%    M = [3 1;2 -4];
%    cPZ_ = M * cPZ;
% 
%    figure; hold on;
%    plot(cPZ,[1,2],'r');
%    plot(cPZ_,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/mtimes, plus

% Author:       Niklas Kochdumper
% Written:      15-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% find the conPolyZono object
[cPZ,M] = findClassArg(factor1,factor2,'conPolyZono');

try

    % normal matrix
    if isnumeric(M)
        
        cPZ.c = M*cPZ.c;

        if ~isempty(cPZ.G)
            cPZ.G = M*cPZ.G;
        end
    
        if ~isempty(cPZ.Grest)
            cPZ.Grest = M*cPZ.Grest;
        end
        
    % use polynomial zonotope method for other cases
    else
        
        pZ = polyZonotope(cPZ.c,cPZ.G,cPZ.Grest,cPZ.expMat);
        
        temp = mtimes(M,pZ);
        
        cPZ = conPolyZono(temp.c,temp.G,temp.expMat,cPZ.A,cPZ.b, ...
                          cPZ.expMat_,temp.Grest,cPZ.id);
        
    end

catch ME
    % note: error has already occured, so the operations below don't have
    % to be efficient

    % already know what's going on...
    if startsWith(ME.identifier,'CORA')
        rethrow(ME);
    end

    % check for empty sets
    if isempty(cPZ)
        return
    end

    % check whether different dimension of ambient space
    equalDimCheck(factor1,factor2);

    % other error...
    rethrow(ME);


end

%------------- END OF CODE --------------