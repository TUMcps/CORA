function [c, Gunred, Gred, indRed] = pickedGeneratorsFast(Z,order)
% pickedGeneratorsFast - Selects generators to be reduced without sorting
%
% Syntax:
%    [c, Gunred, Gred] = pickedGeneratorsFast(Z,order)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    c - center of reduced zonotope
%    Gunred - generators that are not reduced
%    Gred - generators that are reduced
%    indRed - indices that are reduced
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%

% Authors:       Matthias Althoff
% Written:       11-October-2017 
% Last update:   28-October-2017
%                14-March-2019 (vector norm exchanged, remove sort)
%                27-August-2019
%                20-July-2023 (TL, split mink/maxk)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%center
c = center(Z);

%extract generator matrix
G = generators(Z);

%default values
Gunred = [];
Gred = [];
indRed = [];

if ~isempty(G)
    
    %delete zero-length generators
    G = nonzeroFilter(G);

    %number of generators
    [d, nrOfGens] = size(G);
    
    %only reduce if zonotope order is greater than the desired order
    if nrOfGens>d*order

        %number of generators that are not reduced
        nUnreduced = floor(d*(order-1));
        %number of generators that are reduced
        nReduced = nrOfGens - nUnreduced;

        if nReduced == nrOfGens
            % all generators are reduced
            Gred = G;

        else
            %compute metric of generators
            h = vecnorm(G,1,1) - vecnorm(G,Inf,1);

            if nReduced < nUnreduced
                % pick generators with smallest h values to be reduced
                [~,indRed] = mink(h,nReduced);
                idxRed = false(1,nrOfGens);
                idxRed(indRed) = true;
                idxUnred = ~idxRed;

            else
                % pick generators with largest h values to be kept
                [~,indUnred] = maxk(fliplr(h),nUnreduced);
                indUnred = nrOfGens - indUnred + 1; % maintain ordering
                idxUnred = false(1,nrOfGens);
                idxUnred(indUnred) = true;
                idxRed = ~idxUnred;

            end

            % split G accordingly
            Gred = G(:,idxRed);
            Gunred = G(:,idxUnred);
        end

    else
        % no reduction
        Gunred = G;
    end
    
end


% ------------------------------ END OF CODE ------------------------------
