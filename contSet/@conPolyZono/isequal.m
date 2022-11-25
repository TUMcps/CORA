function res = isequal(cPZ1,cPZ2,tol)
% isequal - checks if two conPolyZono are equal
%
% Syntax:  
%    res = isequal(cPZ1,cPZ2)
%    res = isequal(cPZ1,cPZ2,tol)
%
% Inputs:
%    cPZ1 - conPolyZono object
%    cPZ2 - conPolyZono object
%    tol - tolerance (optional)
%
% Outputs:
%    res - boolean whether cPZ1 and cPZ2 are equal
%
% Example: 
%    cPZ1 = conPolyZono([0;0],[1 0 1;0 -1 1],[1 0 2;0 1 1],[1 -0.5], ...
%                       -1,[0 2;1 0],[0.4 0;0.1 1]);
%    cPZ2 = conPolyZono([0;0],[1 1 0;1 0 -1],[2 1 0;1 0 1],[0.5 -1], ...
%                       1,[2 0;0 1],[0 0.4;1 0.1]);
%    isequal(cPZ1,cPZ2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/isequal, zonotope/isequal

% Author:        Niklas Kochdumper
% Written:       27-January-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

    res = false;

    % parse input arguments
    if nargin == 2
        tol = eps;
    end

    % remove redundancies
    cPZ1 = compact(cPZ1);
    cPZ2 = compact(cPZ2);
    
    % compare number of constraints (quick check)
    if size(cPZ1.A,1) ~= size(cPZ2.A)
        return;
    end
    
    % compare number of generators (quick check)
    if size(cPZ1.G,2) ~= size(cPZ2.G,2) || ... 
       size(cPZ1.A,2) ~= size(cPZ2.A,2) || ...
       size(cPZ1.Grest,2) ~= size(cPZ2.Grest,2)
       return 
    end
    
    % compare identifier vectors
    temp1 = sort(cPZ1.id); temp2 = sort(unique([cPZ1.id;cPZ2.id]));
    E1 = cPZ1.expMat; E2 = cPZ2.expMat; 
    E1_ = cPZ1.expMat_; E2_ = cPZ2.expMat_;
    
    if length(temp1) ~= length(temp2) || ~all(temp1 == temp2)
       return;
    elseif ~all(cPZ1.id == cPZ2.id)
       [~,E1,E2] = mergeExpMatrix(cPZ1.id,cPZ2.id,cPZ1.expMat,cPZ2.expMat);
       if ~isempty(cPZ1.A)
        [~,E1_,E2_] = mergeExpMatrix(cPZ1.id,cPZ2.id,cPZ1.expMat_, ...
                                     cPZ2.expMat_);
       end
    end
    
    % compare exponent matrices
    if ~(all(all(abs(E1 - E2) < tol)))
        return
    end

    % compare dependent generators
    if ~(all(all(abs(cPZ1.G - cPZ2.G) < tol)))
        return
    end
    
    % compare constraints
    if ~isempty(cPZ1.A)
       
        % compare constraint exponent matrices
        if ~(all(all(abs(E1_ - E2_) < tol)))
            return
        end
        
        % compare constraints (consider that f(x) = 0 is identical to -f(x)
        % = 0)
        A1 = [cPZ1.b,cPZ1.A]; A2 = [cPZ2.b,cPZ2.A];
        
        for i = 1:size(A1,1)
           if ~(all(all(abs(A1(i,:) - A2(i,:)) < tol)))
              if ~(all(all(abs(-A1(i,:) - A2(i,:)) < tol)))
                  return;
              end
           end
        end
    end

    % compare center and independent generators
    Z1 = zonotope(cPZ1.c,cPZ1.Grest);
    Z2 = zonotope(cPZ2.c,cPZ2.Grest);
    res = isequal(Z1,Z2);
end

%------------- END OF CODE --------------