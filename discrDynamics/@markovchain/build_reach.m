function [Obj,tp]=build_reach(Obj,R,iInput,iState,tranFrac)
% build - Builds the transition matrices of the Markov chains using 
% reachability analysis.
%
% Syntax:  
%    [Obj,tp] = build(Obj,R,iInput,iState,tranFrac)
%
% Inputs:
%    Obj - Markov chain object
%    R - reachSet object
%    iInput - discrete input
%    iState - initial discrete state
%    tranFrac - 
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
%               24-July-2020
% Last revision:---

%------------- BEGIN CODE --------------

%get probabilities from actual segment to all segments
tp = transitionProbability_reach(R,tranFrac.TP,Obj.field);

%load transition probability from actual segment to reachable segments in
%Transition Matrix T
actualSegmentNr = iState;
Obj.T.T{iInput}(:,actualSegmentNr+1)=sparse(tp.T);
Obj.T.OT{iInput}(:,actualSegmentNr+1)=sparse(tp.OT);

%outside probabilities should stay outside
Obj.T.T{iInput}(1,1)=1;
Obj.T.OT{iInput}(1,1)=1;

%------------- END OF CODE --------------