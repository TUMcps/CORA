function [lcEvolProb]=laneChangeEvolution(pChange,lcEvolProb,markovChainSpec)
% laneChangeEvolution - updates the probability distribution for the state
% of the lane change from the beginning up to the end of the lane change
% process. The probabilities for driving on the left and right lane is also
% implicitly specified.
%
% Syntax:  
%    [lcEvolProb]=laneChangeEvolution(pChange,lcEvolProb)
%
% Inputs:
%    pChange - probability that a lane change has started
%    lcEvolProb - matrix which specifies the status of the lane changes
%    over time
%    markovChainSpec - struct specifying Markov chain properties
%
% Outputs:
%    lcEvolProb - matrix which specifies the status of the lane changes
%    over time
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      03-December-2008 
% Last update:  16-October-2009
% Last revision:---

%------------- BEGIN CODE --------------

%set duration for a lane change
lcDuration=5; %sec
%intermediate states of a lane change 
lcSteps=ceil(lcDuration/markovChainSpec.timeStep);

if isempty(lcEvolProb)
    lcEvolProb=zeros(1,lcSteps+1);
    
    %insert actual lane change probability
    lcEvolProb(1,1)=1-pChange; %prob of driving left
    lcEvolProb(1,2)=pChange; %prob of starting lane change    
else
    %obtain time step
    iStep=length(lcEvolProb(:,1));

    %insert actual lane change probability
    lcEvolProb(iStep+1,1)=lcEvolProb(iStep,1)*(1-pChange); %prob of driving left
    lcEvolProb(iStep+1,2)=lcEvolProb(iStep,1)*pChange; %prob of starting lane change

    %forward probabilities to next lc phase
    for i=2:lcSteps
        lcEvolProb(iStep+1,i+1)=lcEvolProb(iStep,i);
    end

    %add probability of driving on right lane
    lcEvolProb(iStep+1,i+1)=lcEvolProb(iStep+1,i+1)+lcEvolProb(iStep,i+1);
end


%------------- END OF CODE --------------