function [Rnext,IH] = postReach(obj,Rinit,R_tp,c)
% postReach - computes the reachable continuous set for the first time
% interval as a postprocessing step
%
% Syntax:  
%    [Rnext] = postReach(obj,Rinit,R_tp)
%
% Inputs:
%    obj - linearSys object
%    Rinit - initial reachable set
%    R_tp - reachable set at the next point in time
%
% Outputs:
%    Rnext - next reachable set 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      03-May-2011
% Last update:  25-July-2016 (intervalhull replaced by interval)---
% Last revision:---

%------------- BEGIN CODE --------------


%load data from object structure
F=obj.F;
inputCorr=0; %add input correction later
RV=obj.RV;


%time interval solution
R_err = F*(Rinit+(-c))+inputCorr;
R_ti=enclose(Rinit,R_tp)+R_err;

%write results to reachable set struct Rfirst
Rnext.tp=R_tp+RV;
Rnext.ti=R_ti+RV;

%compute enclosing hull
IH_init = interval(Rinit);
IH_tp = interval(R_tp);
IH_err = interval(R_err+RV);

IH = enclose(IH_init,IH_tp) + IH_err;

%IH_alt = interval(R_ti);


%------------- END OF CODE --------------