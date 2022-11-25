function [ExpMatNew,Gnew] = removeRedundantExponents(ExpMat,G)
% removeRedundantExponents - add up all generators that belong to terms
%                            with identical exponents 
%
% Syntax:  
%    [ExpMatNew, Gnew] = removeRedundantExponents(ExpMat, G)
%
% Inputs:
%    ExpMat - matrix containing the exponent vectors
%    G - generator matrix
%
% Outputs:
%    ExpMat - modified exponent matrix
%    G - modified generator matrix
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:        Niklas Kochdumper
% Written:       25-June-2018 
% Last update:   21-April-2020 (remove zero-length generators)
% Last revision: ---

%------------- BEGIN CODE --------------

    % remove zero-length generators
    idxD = any(G,1);
    
    % skip if all non-zero
    if ~all(idxD)
        if all(~idxD)
            ExpMatNew = zeros(size(ExpMat,1),1);
            Gnew = zeros(size(G,1),1);
            return;
        else
            G = G(:,idxD);
            ExpMat = ExpMat(:,idxD);
        end
    end

    % add hash value of the exponent vectors to the exponent matrix 
    temp = 1:1:size(ExpMat,1);
    rankMat = [(temp*ExpMat)', ExpMat'];
    
    % sort the exponent vectors according to the hash value
    [~, ind] = sortrows(rankMat);
    ExpMatTemp = ExpMat(:,ind);
    Gtemp = G(:,ind);
    
    % initialization
    counterNew = 1;
    ExpMatNew = zeros(size(ExpMat));
    Gnew = zeros(size(G));
    
    % first entry
    ExpMatNew(:,counterNew) = ExpMatTemp(:,1);
    Gnew(:,counterNew) = Gtemp(:,1);
    
    % loop over all exponent vectors
    for counter = 2:size(ExpMatTemp,2)
        
        if all(ExpMatNew(:,counterNew) == ExpMatTemp(:,counter))           
            Gnew(:,counterNew) = Gnew(:,counterNew) + Gtemp(:,counter);
        else
            counterNew = counterNew + 1;
            Gnew(:,counterNew) = Gtemp(:,counter);
            ExpMatNew(:,counterNew) = ExpMatTemp(:,counter);
        end
        
    end 
    
    % truncate exponent and generator matrix
    ExpMatNew = ExpMatNew(:,1:counterNew);
    Gnew = Gnew(:,1:counterNew);

end

%------------- END OF CODE --------------