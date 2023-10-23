function [Gred]=generatorVolumeFilter(G,rem)
% generatorVolumeFilter - filters out generators by finding the
% combinations returning the biggest volume
%
% Syntax:
%    [Gred]=generatorVolumeFilter(G,rem)
%
% Inputs:
%    G - matrix of generators
%    rem - number of remaining generators
%
% Outputs:
%    Gred - cell array of generators
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Matthias Althoff
% Written:       12-September-2008
% Last update:   19-July-2010
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%determine generators by volume maximation:
%possible combinations of n=dim generators from all generators
[rows,cols]=size(G);
comb = combinator(cols,rows,'c');
nrOfComb=length(comb(:,1));

for i=1:nrOfComb
    try
        %obtain Gpicked
        Gpicked=G(:,comb(i,:));
        parallelogramVol(i)=abs(det(Gpicked));
    %     %check rank of picked generators
    %     if rank(Gpicked)<rows
    %         Gpicked=Gpicked+1e-6*eye(rows);
    %     end
    
    catch
        parallelogramVol(i)=0;
        disp('parallelogram volume could not be computed');
    end
end

% obtain indices corresponding to the largest values
[~,index]=maxk(vol, rem);

%store the generator combinations in cells
for i=1:length(index)
    generatorIndices=comb(index(i),:);
    Gred{i}=G(:,generatorIndices);
end


% ------------------------------ END OF CODE ------------------------------
