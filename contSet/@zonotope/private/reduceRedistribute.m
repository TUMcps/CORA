function [Zred]=reduceRedistribute(Z,order)
% reduceRedistribute - Reduce remaining generators of a zonotope
% so that its order stays below a specified limit 
%
% Syntax:
%    [Zred]=reduceRedistribute(Z,order)
%
% Inputs:
%    Z - zonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    Z - zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       07-September-2012 
% Last update:   16-March-2019 (vnorm replaced, sort removed)
%                27-August-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%initialize Z_red
Zred=Z;

%extract center and generator matrix
c = center(Z);
G = generators(Z);

%dimension and number of generators
[d, nrOfGens] = size(G);

if ~isempty(G)

    %only reduce if zonotope order is greater than the desired order
    if nrOfGens>d*order

        %compute metric of generators (shortest generators)
        h=vecnorm(G);
        
        %remove elements of length less than 1e-10
        [~, fInd] = find(h<max(max(G))*1e-6);
        G(:,fInd) = [];
        
        %only reduce if zonotope order is greater than the desired order
        if length(G(1,:))>d*order
            
            %compute metric of generators (shortest generators)
            h=vecnorm(G);

            %number of generators that are not reduced
            nUnreduced=floor(d*order);
            %number of generators that are reduced
            nReduced=length(G(1,:))-nUnreduced;

            %pick generators that are reduced
            [~,ind] = mink(h,nReduced);
            pickedGens=G(:,ind);
            
            %unreduced generators
            indRemain = setdiff(1:nrOfGens, ind);
            Gunred=G(:,indRemain);
            
            %scale generators in G for compensation
            Gnew = aux_generatorScaling(Gunred, pickedGens);
            %Gold = aux_generatorScaling_old(Gunred, pickedGens);

            %build reduced zonotope
            Zred.c = c;
            Zred.G = Gnew;
        end

    end
end


% Auxiliary functions -----------------------------------------------------

function Gnew = aux_generatorScaling(Grem, Gdel)

%dim 
d = length(Grem(:,1));

%remove too small generators
scaleFactor = vecnorm(Grem);
[~,ind] = find(scaleFactor>0);

%normalize remaining generators
for i=1:length(Grem(:,1))
    Gnorm(i,:) = Grem(i,ind)./scaleFactor(ind);
end

%check alignment of each generator in Gdel
scale = ones(length(Grem(1,:)),1);

%get frame out of most and least aligned generators
perpendicularInd_pre = aux_pickPerpendicular(Gnorm,d);


for i=1:length(Gdel(1,:))
    prod = abs(Gdel(:,i)'*Gnorm);
    
    %get frame out of most and least aligned generators
    %remove largest value
    [~,pickedInd]=max(prod);
    perpendicularInd = setdiff(perpendicularInd_pre,pickedInd);
    if length(perpendicularInd) == d
        [~,ind]=max(abs(Gnorm(:,pickedInd)'*Gnorm(:,perpendicularInd)));
        perpendicularInd(ind) = [];
    end
    chosenInd = [perpendicularInd,pickedInd];
    frame = Gnorm(:,chosenInd);
    scaleFactorSort(:,1) = scaleFactor(chosenInd);
    
    %add to scaling
    %addedScaling_old = abs(pinv(frame)*Gdel(:,i))./scaleFactorSort;
    addedScaling = abs(frame\Gdel(:,i))./scaleFactorSort;
    if any(addedScaling>1e10)
        disp('stop!!!');
    end
%     for i=1:length(Grem(:,1))
%         addedG(i,:) = Grem(i,chosenInd).*addedScaling';
%     end
    scale(chosenInd) = scale(chosenInd) + addedScaling;
end

%scale remaining generators
for i=1:length(Grem(:,1))
    Gnew(i,:) = Grem(i,:).*scale';
end


%pick n-1 perpendicular generators
function perpendicularInd = aux_pickPerpendicular(Gnorm,dim)

%which generatpors are not least aligned with all other generators?
alignmentMat = abs(Gnorm'*Gnorm);

%remove diagonals
alignmentMat = alignmentMat - diag(diag(alignmentMat));

%least maximum entry?
finalInd = 1:length(alignmentMat);
[elements,indices] = max(alignmentMat);
[~,maxInd] = max(elements);
remInd1 = indices(maxInd);
remInd2 = finalInd(remInd1);

%remove rows and columns
while length(finalInd)>dim
    %remove elements form alignment matrix
    alignmentMat(:,remInd1) = [];
    alignmentMat(remInd1,:) = [];
    finalInd = setdiff(finalInd, remInd2);
    
    %new check
    [elements,indices] = max(alignmentMat);
    [~,maxInd] = max(elements);
    remInd1 = indices(maxInd);
    remInd2 = finalInd(remInd1);
end

%choose smallest values
perpendicularInd = finalInd;

% %pick n-1 perpendicular generators
% function perpendicularInd_old = aux_pickPerpendicular(Gnorm,pickedInd,dim)
% 
% %which generatpors are not least aligned with all other generators?
% alignmentMat = abs(Gnorm'*Gnorm);
% 
% %remove diagonals
% alignmentMat = alignmentMat - diag(diag(alignmentMat));
% 
% %least maximum entry?
% [elements,indices] = sort(max(alignmentMat));
% 
% %remove picked Ind
% indRem = find(indices == pickedInd);
% indices(indRem) = [];
% 
% %choose smallest values
% perpendicularInd = indices(1:(dim-1));


function Gnew = aux_generatorScaling_old(Grem, Gdel)

%dim 
d = length(Grem(:,1));

%remove too small generators
scaleFactor = vecnorm(Grem);
[~,ind] = find(scaleFactor>0);

%normalize remaining generators
for i=1:length(Grem(:,1))
    Gnorm(i,:) = Grem(i,ind)./scaleFactor(ind);
end

%check alignment of each generator in Gdel
scale = ones(length(Grem(1,:)),1);
for i=1:length(Gdel(1,:))
    prod = abs(Gdel(:,i)'*Gnorm);
    [~,indices]=sort(prod);
    
    %get frame out of most and least aligned generators
    chosenInd = [indices(1:d-1),indices(end)];
    frame = Grem(:,chosenInd);
    
    %add to scaling
    addedScaling = abs(pinv(frame)*Gdel(:,i));
%     if any(addedScaling>1.5)
%         disp('stop!!!');
%     end
    scale(chosenInd) = scale(chosenInd) + addedScaling;
end

%scale remaining generators
for i=1:length(Grem(:,1))
    Gnew(i,:) = Grem(i,:).*scale';
end


% ------------------------------ END OF CODE ------------------------------
