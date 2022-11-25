function [Rfirst,options] = initReach(obj,Rinit,options)
% initReach - computes the reachable continuous set for the first time step
%
% Syntax:  
%    [obj,Rfirst] = initReach(obj,Rinit,options)
%
% Inputs:
%    obj - linProbSys object
%    Rinit - initial reachable set struct
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linIntSys object
%    Rfirst - first reachable set struct
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      13-September-2007 
% Last update:  06-October-2007
%               08-September-2009
%               25-July-2016 (intervalhull replaced by interval)
%               17-July-2020 (return values changed)
% Last revision:---

%------------- BEGIN CODE --------------


% compute exponential matrix
obj = pexpm(obj,options);
% compute time interval error (tie)
obj = tie(obj,options);
% compute reachable set due to input
obj = inputSolution(obj,options);
%change the time step
obj.taylor.timeStep=options.timeStep;

%load data from object structure
eAt=obj.taylor.eAt;
F=obj.taylor.F;
Rinput=obj.taylor.Rinput;
Rtrans=obj.taylor.Rtrans;
inputCorr=obj.taylor.inputCorr;
pRinput=obj.taylor.pRinput;

%first time step homogeneous solution
Rhom_tp=eAt*Rinit;

%compute auxiliary enclosing probabilistic hull
Raux=enclose(Rinit,obj.A*options.timeStep,Rtrans);
%enclarge due to time interval 
Rhom=Raux+F*zonotope(Rinit)+inputCorr+(-Rtrans);

%particulate solution
Rpar=Rinput+pRinput;

%total solution
Rtotal=Rhom+zonotope(Rpar); %<-- mSigma-bound for first particulate reachable set!
Rtotal_tp=Rhom_tp+zonotope(Rpar);

%save homogeneous and particulate solution
options.Rhom=Rhom;
%options.Rhom=Rhom_tp; %<--change for time point solution
options.Raux_det=zonotope(Rinput,0);
options.Rpar_det=interval(options.Raux_det);
options.Rpar_prob=zonotope(pRinput);
%options.Rpar_prob=pRinput; %<--change for time point solution

%write results to reachable set struct Rfirst
Rfirst.tp=Rtotal_tp;
Rfirst.ti=Rtotal;

%------------- END OF CODE --------------