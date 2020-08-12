function [id,expMat1,expMat2] = mergeExpMatrix(id1,id2,expMat1,expMat2)
% mergeExpMatrix - Merge the ID-vectors of two polyZonotope objects
%                  and adapte the exponent matrices accordingly
%
% Syntax:  
%    [id,expMat1,expMat2] = mergeExpMatrix(id1,id2,expMat1,expMat2)
%
% Inputs:
%    id1 - ID-vector of the first polynomial zonotope
%    id2 - ID-vector of the second polynomial zonotope
%    expMat1 - exponent matrix of the first polynomial zonotope
%    expMat2 - exponent matrix of the second polynomial zonotope
%
% Outputs:
%    id - merged ID-vector
%    expMat1 - adapted exponent matrix of the first polynomial zonotope
%    expMat2 - adapted exponent matrix of the second polynomial zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      25-June-2018 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    L1 = length(id1);
    L2 = length(id2);

    % ID vectors are identical
    if L1 == L2 && all(id1 == id2)
       
        id = id1;
        
    % ID vectors not identical -> MERGE
    else
        
        % merge the two sets
        id = id1;
        ind2 = zeros(size(id2));
        for i = 1:length(id2)
           ind = find(id == id2(i));
           if isempty(ind)
              id = [id;id2(i)];
              ind2(i) = length(id);
           else
              ind2(i) = ind;
           end
        end
        
        % construct the new exponent matrices
        L = length(id);
        
        expMat1 = [expMat1;zeros(L-L1,size(expMat1,2))];
        
        temp = zeros(L,size(expMat2,2));
        temp(ind2,:) = expMat2;
        expMat2 = temp;
    end
end

%------------- END OF CODE --------------