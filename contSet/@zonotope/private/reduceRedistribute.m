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

% Author:       Matthias Althoff
% Written:      07-September-2012 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%initialize Z_red
Zred=Z;

%get Z-matrix from zonotope Z
Zmatrix=get(Z,'Z');

%extract generator matrix
G=Zmatrix(:,2:end);

if ~isempty(G)

    %determine dimension of zonotope
    dim=length(G(:,1));

    %only reduce if zonotope order is greater than the desired order
    if length(G(1,:))>dim*order

        %compute metric of generators (shortest generators)
        h=vnorm(G,1,2);
        
        %remove elements of length less than 1e-10
        [fElem, fInd] = find(h<max(max(G))*1e-6);
        G(:,fInd) = [];
        
        %only reduce if zonotope order is greater than the desired order
        if length(G(1,:))>dim*order
            
            %compute metric of generators (shortest generators)
            h=vnorm(G,1,2);
            [elements,indices]=sort(h);

            %number of generators that are not reduced
            nUnreduced=floor(dim*(order));
            %number of generators that are reduced
            nReduced=length(G(1,:))-nUnreduced;
            
            %unreduced generators
            Gunred=G(:,indices((nReduced+1):end));

            %pick generators that are reduced
            pickedGenerators=G(:,indices(1:nReduced));
            %scale generators in G for compensation
            Gnew = generatorScaling(Gunred, pickedGenerators);
            %Gold = generatorScaling_old(Gunred, pickedGenerators);

            %build reduced zonotope
            Zred.Z=[Zmatrix(:,1),Gnew];
        end

    end
end


function Gnew = generatorScaling(Grem, Gdel)

%dim 
dim = length(Grem(:,1));

%remove too small generators
scaleFactor = vnorm(Grem,1,2);
[val,ind] = find(scaleFactor>0);

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
[maxElem,maxInd] = max(elements);
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
    [maxElem,maxInd] = max(elements);
    remInd1 = indices(maxInd);
    remInd2 = finalInd(remInd1);
end

%choose smallest values
perpendicularInd = finalInd;

% %pick n-1 perpendicular generators
% function perpendicularInd_old = pickPerpendicular(Gnorm,pickedInd,dim)
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


function Gnew = generatorScaling_old(Grem, Gdel)

%dim 
dim = length(Grem(:,1));

%remove too small generators
scaleFactor = vnorm(Grem,1,2);
[val,ind] = find(scaleFactor>0);

%normalize remaining generators
for i=1:length(Grem(:,1))
    Gnorm(i,:) = Grem(i,ind)./scaleFactor(ind);
end

%check alignment of each generator in Gdel
scale = ones(length(Grem(1,:)),1);
for i=1:length(Gdel(1,:))
    prod = abs(Gdel(:,i)'*Gnorm);
    [elements,indices]=sort(prod);
    
    %get frame out of most and least aligned generators
    chosenInd = [indices(1:dim-1),indices(end)];
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


%------------- END OF CODE --------------