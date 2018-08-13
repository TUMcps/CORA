function matV = dominantVertices(matI, maxNumber)
% dominantVertices - computes the dominant vertices of an interval matrix
%
% Syntax:  
%    matV = dominantVertices(matI, maxNumber)
%
% Inputs:
%    matI - interval matrix
%    maxNumber - maximum number of dominant vertices
%
% Outputs:
%    matV - cell array of matrix vertices
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      02-July-2010 
% Last update:  19-July-2010
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%conversion to an interval
IH = interval(matI);

%get edge lengths
edgeLen = 2*rad(IH);

%sort edge lengths
[edgeLenSorted, ind] = sort(edgeLen,'descend');

%set default value
if maxNumber==inf
    unlimited=1;
    nrOfParameters = length(find(edgeLen>0));
    maxNumber = 2^nrOfParameters;
else
    unlimited=0;
end

%start at center 
V = center(IH);

%add new vertices
counter = 1;
nrOfVertices=1;
while ((nrOfVertices*2)<(maxNumber+2)) && (counter<=length(edgeLen))
    %current index
    currInd = ind(counter);
    %obtain Vleft vertices
    Vleft = V;
    Vleft(currInd,:) = V(currInd,:) - 0.5*edgeLenSorted(counter);
    %obtain Vright vertices
    Vright = V;
    Vright(currInd,:) = V(currInd,:) + 0.5*edgeLenSorted(counter);
    %combine vertices
    V = [Vleft,Vright];
    %update counter, nrOfVertices
    counter = counter+1;
    nrOfVertices = nrOfVertices*2;
end

%convert vertices to matrix vertices
if unlimited
    matV=cell(length(V(1,:)),1);
else
    matV=cell(length(V(1,:))+2,1);
end
for i=1:length(V(1,:))
    matV{i}=vec2mat(V(:,i));
end

%sup and inf matrices are also dominant!
if ~unlimited
    matV{i+1} = infimum(matI.int);
    matV{i+2} = supremum(matI.int);
end

%------------- END OF CODE --------------