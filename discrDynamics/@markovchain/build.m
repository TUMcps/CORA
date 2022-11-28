function [Obj,tp] = build(Obj,finalStateMat,iInput,stateField)
% build - Builds the transition matrices of the Markov chains using 
% simulation.
%
% Syntax:  
%    [Obj,tp] = build(Obj,finalStateMat,iInput,stateField)
%
% Inputs:
%    Obj - Markov chain object
%    finalStateMat -
%    iInput - discrete input
%    stateField - partition of the continuous state space
%
% Outputs:
%    Obj - Markov chain object
%    tp - transition probability struct
%
% Example: 
%    -
%
% Other m-files required: interval (constructor)
% Subfunctions: none
% MAT-files required: none
%
% See also: vertices, polytope

% Author:       Matthias Althoff
% Written:      15-September-2006
% Last update:  28-September-2006
%               16-August-2007
%               23-November-2007
%               21-April-2009
%               16-June-2009
% Last revision:---

%------------- BEGIN CODE --------------

%compute transition probabilities
[tp.T]=transitionProbability(finalStateMat.T,Obj.field);
[tp.OT]=transitionProbability(finalStateMat.OT,Obj.field);

%load transition probability from actual segment to reachable segments in
%Transition Matrix T
%actualSegmentNr=get(Obj.field,'actualSegmentNr');
actualSegmentNr = stateField; %I THINK -- AP
Obj.T.T{iInput}(:,actualSegmentNr+1)=sparse(tp.T);
Obj.T.OT{iInput}(:,actualSegmentNr+1)=sparse(tp.OT);

%outside probabilities should stay outside
Obj.T.T{iInput}(1,1)=1;
Obj.T.OT{iInput}(1,1)=1;

%------------- END OF CODE --------------