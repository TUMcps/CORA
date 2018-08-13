function [qZ]=reduceRedistribute(qZ,order)
% reduceRedistribute - Reduce remaining generators of a quadratic zonotope
% so that its order stays below a specified limit 
%
% Syntax:  
%    [qZred]=reduceRedistribute(qZ,order)
%
% Inputs:
%    qZ - quadZonotope object
%    order - desired order of the zonotope
%
% Outputs:
%    qZ - reduced quadZonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      06-September-2012 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


if ~isempty(qZ.Grest)

    %determine dimension of zonotope
    dim=length(qZ.Grest(:,1));


    %only reduce if zonotope order is greater than the desired order
    if length(qZ.Grest(1,:))>dim*order

        %compute metric of generators (shortest generators)
        h=vnorm(qZ.Grest,1,2);
        
        %remove elements of length less than 1e-10
        [~, fInd] = find(h<max(max(qZ.Grest))*1e-6);
        qZ.Grest(:,fInd) = [];
        
        %only reduce if zonotope order is greater than the desired order
        if length(qZ.Grest(1,:))>dim*order
            
            %compute metric of generators (shortest generators)
            h=vnorm(qZ.Grest,1,2);
            [~,indices]=sort(h);

            %number of generators that are not reduced
            nUnreduced=floor(dim*order);
            %number of generators that are reduced
            nReduced=length(qZ.Grest(1,:))-nUnreduced;

            %pick generators that are reduced
            pickedGenerators=qZ.Grest(:,indices(1:nReduced));
            %scale generators in G for compensation
            Gnew = generatorScaling(qZ.G, pickedGenerators);
            qZ.G = Gnew;

            %unreduced generators
            qZ.Grest=qZ.Grest(:,indices((nReduced+1):end));
        end

    end
end


function Gnew = generatorScaling(Grem, Gdel)

%dim 
dim = length(Grem(:,1));

%remove too small generators
scaleFactor = vnorm(Grem,1,2);
[~,ind] = find(scaleFactor>0);

%normalize remaining generators
for i=1:length(Grem(:,1))
    Gnorm(i,:) = Grem(i,ind)./scaleFactor(ind);
end

%check alignment of each generator in Gdel
scale = ones(length(Grem(1,:)),1);

%get frame out of most and least aligned generators
perpendicularInd_pre = pickPerpendicular(Gnorm,dim);


for i=1:length(Gdel(1,:))
    prod = abs(Gdel(:,i)'*Gnorm);
    [~,indices]=sort(prod);
    
    %get frame out of most and least aligned generators
    %remove picked Ind
    pickedInd = indices(end);
    perpendicularInd = setdiff(perpendicularInd_pre,pickedInd);
    if length(perpendicularInd) == dim
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
function perpendicularInd = pickPerpendicular(Gnorm,dim)

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



%------------- END OF CODE --------------