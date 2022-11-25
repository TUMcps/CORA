function [T] = roadTransitionCopying(tp,T,field)
% roadTransitionCopying - Copies transitions describing motion of traffic
% participants on a road; the transitions can be copied as the transition
% probabilities are independent of the position of the traffic participant
% on the road; this function is only for road objects discretized in
% position and velocity
%
% Syntax:  
%    [T] = roadTransitionCopying(tp,T,field)
%
% Inputs:
%    tp - transition probability vector for a specific input and velocity
%    interval
%    T - Transition matrix of the Markov chain
%    field - partition object of the discretized state space
%
% Outputs:
%    T - Transition matrix of the Markov chain
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written: 27-March-2008
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

%generate auxiliary variable for the transition probabilities
trans=tp;

%delete outside transition probabilities
trans(1)=[];

%find nonzero transition probabilities
cellNr=find(trans);

%get nr of segments
%nrOfSegments=get(field,'nrOfSegments'); %<--AP
nrOfSegments = field.nrOfSegments;

if ~isempty(cellNr)

    %obtain subscript matrix
    subscriptMatrix=cellSegments(field,cellNr');  
    
    %determine maximum reached position segment
    maxPos=max(subscriptMatrix(:,1)); % <-- possible problem, AP

    %number of direct copies
    dCopies=nrOfSegments(1)-maxPos;

    %copy transitions
    tp=sparse(tp);
    indices=find(tp(2:end))+1;
    for i=1:dCopies
        T(1,end+1)=tp(1);
        T(i+indices,end)=sparse(tp(indices));
    end

    %copy transitions that may hit outside area
    for iPos=1:(maxPos-1)
        %transition probability to outside
        T(1,end+1)=tp(1);
        for iInd=1:length(indices)
            iCell=indices(iInd)-1;
            %obtain subscript vector
            subscriptVector=cellSegments(field,iCell); 
            if subscriptVector(1)<=maxPos-iPos
                %copy transition probabilities
                T(dCopies+iPos+iCell+1,end)=sparse(tp(iCell+1));
            else
                %increase outside probability
                T(1,end)=T(1,end)+sparse(tp(iCell+1));
            end
        end
    end
else
    T(nrOfSegments+1,nrOfSegments+1)=0;
end

%------------- END OF CODE --------------