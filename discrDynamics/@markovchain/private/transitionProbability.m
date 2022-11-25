function tp = transitionProbability(finalStateMat,field)
% transitionProbability - Calculate the transition probability from the 
% actual cell to the reachable cells.
%
% Syntax:  
%    tp =transitionProbability(finalStateMat,field)
%
% Inputs:
%    finalStateMat - 
%    field - partition of the state space 
%
% Outputs:
%    tp - transition probability struct
%
% Example: 
%    -
%
% Other m-files required: tbd
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      15-September-2006
% Last update:  09-October-2006
%               26-March-2008
%               16-June-2009
%               14-August-2018
% Last revision:---

%------------- BEGIN CODE --------------


%initialize--------------------------------------------------------
nrOfStates = nrOfCells(field);
tp(1:(nrOfStates+1),1)=0; %tp: transition probability
%------------------------------------------------------------------

%get total number of final states----------------------------------
nrOfFinalStates=length(finalStateMat(:,1));
%------------------------------------------------------------------

%get cell indices of final states----------------------------------
finalCell = zeros(1,nrOfFinalStates);
for iPoint=1:nrOfFinalStates
    aux=intersectingCells(field,finalStateMat(iPoint,:));
    finalCell(iPoint)=aux(1); %if at border of cells, choose first cell
end
%------------------------------------------------------------------

%calculate probabilities for each cell by: number of trajectories entering
%the goal state devided by the number of all trajectories
while ~isempty(finalCell)
    indices=find(finalCell == finalCell(1));
    tp(finalCell(1)+1,1)=length(indices)/nrOfFinalStates;
    finalCell(indices)=[];
end
%------------------------------------------------------------------

%------------- END OF CODE --------------