function [obj,stateField]=build4road_reach(obj,R,inputNr,stateField)
% build4road - Specialized function to build Markov chains, based on the
% fact that reachable sets for road vehicles are independent of their
% position
%
% Syntax:  
%    [obj,stateField]=build4road(MC,HA,inputNr,stateField)
%
% Inputs:
%    obj - Markov chain object
%    HA - Hybrid automaton
%    inputNr - input Nr of the current discrete input
%    stateField - partition object of the discreteized state space
%
% Outputs:
%    obj - Markov chain object
%    stateField - partition object of the discreteized state space
%
% Example: 
%   ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      27-March-2007 
% Last update:  21-April-2009
%               24-July-2020
% Last revision: ---

%------------- BEGIN CODE --------------

%compute transition fraction for the time point and the time solution
tranFrac=transitionFraction(R);

%build Markov chain for the current initial cell
[obj,tp]=build_reach(obj,R,inputNr,stateField,tranFrac);

%speed up Markov chain generation by copying transition probabilities
[obj.T.T{inputNr}] = roadTransitionCopying(tp.T,obj.T.T{inputNr},obj.field);
[obj.T.OT{inputNr}] = roadTransitionCopying(tp.OT,obj.T.OT{inputNr},obj.field);

%------------- END OF CODE --------------