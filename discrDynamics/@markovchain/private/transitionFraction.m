function [tranFrac] = transitionFraction(R,Rcont,TP)
% transitionFraction - Computes by which probability a transition has
% occured from deceleration to standstill or from acceleration to speed
% limit. The probability is computed straightforward for the time point
% solution; the time interval solution is computed in an approximate manner
% (see AT paper draft)
%
% Syntax:  
%    [tranFrac] = transitionFraction(R,Rcont,tmin,tmax,tend)
%
% Inputs:
%    R - reachable set
%    Rcont - reachable set that has not been intersected with the guard set
%    TP - time struct for initial, last time point in invariant as well as 
%    first and last time point of guard intersection
%
% Outputs:
%    tranFrac - struct for transition probabilities
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
% Written:      21-April-2009
% Last update:  31-July-2017
% Last revision: ---

%------------- BEGIN CODE --------------

%check if a transition has occured
if length(R.T)==1
    tranFrac.TP{1}=1;
    tranFrac.TI{1}=1;    
else

    %determine the transition fraction for the time point solution

    %compute volume of intersected reachable set
    volInt=volume(R.T{1}{end});
    %compute volume of reachable set without intersection
    volOrig=volume(polytope(Rcont.T{1}{end}));   %use polytope conversion as this has been applied to R, too.
    %compute ratio
    ratio=volInt/volOrig;
    
    %set transition fraction for time points
    tranFrac.TP{1}=ratio;
    tranFrac.TP{2}=1-ratio;
    
    %approximate time interval ratio
    ratioTI=TP{1}.tMax/TP{2}.tMax;
    
    %set transition fraction for time intervals
    tranFrac.TI{1}=ratioTI;
    tranFrac.TI{2}=1-ratioTI;    
end


%------------- END OF CODE --------------