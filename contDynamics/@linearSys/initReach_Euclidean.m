function [Rfirst,options] = initReach_Euclidean(obj,Rinit,options)
% reach - computes the reachable continuous set for the first time step in
% the untransformed space
%
% Syntax:  
%    [Rfirst,options] = initReach_Euclidean(obj,Rinit,options)
%
% Inputs:
%    obj - linearSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rfirst - first reachable set 
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author: Matthias Althoff
% Written:      07-May-2007 
% Last update:  03-January-2008
%               04-May-2009
%               29-June-2009
%               08-August-2011
%               25-July-2016 (intervalhull replaced by interval)
%               06-April-2017
%               28-October-2017
%               07-November-2018
% Last revision: ---

%------------- BEGIN CODE --------------

% compute exponential matrix
obj = exponential(obj,options);
% compute time interval error (tie)
obj = tie(obj,options);
% compute reachable set due to input
obj = inputSolution(obj,options);
%change the time step
obj.taylor.timeStep=options.timeStep;

%compute reachable set of first time interval
eAt=expm(obj.A*options.timeStep);
%save data to object structure
obj.taylor.eAt=eAt;

F=obj.taylor.F;
RV=obj.taylor.RV;
inputCorr=obj.taylor.inputCorr;
if iscell(obj.taylor.Rtrans)
    Rtrans=obj.taylor.Rtrans{1};
else
    Rtrans=obj.taylor.Rtrans;
end


%first time step homogeneous solution
Rhom_tp=eAt*Rinit + Rtrans;
if isa(Rinit,'polyZonotope')
    Rhom=enclose(Rinit,Rhom_tp)+F*zonotope(Rinit)+inputCorr;
elseif isa(Rinit,'zonoBundle') 
    Rhom=enclose(Rinit,Rhom_tp)+F*Rinit.Z{1}+inputCorr;
else
    Rhom=enclose(Rinit,Rhom_tp)+F*Rinit+inputCorr;
end

%reduce zonotope
Rhom=reduce(Rhom,options.reductionTechnique,options.zonotopeOrder);
Rhom_tp=reduce(Rhom_tp,options.reductionTechnique,options.zonotopeOrder);
RV=reduce(RV,options.reductionTechnique,options.zonotopeOrder);

%save homogeneous and particulate solution
options.Rhom=Rhom;
options.Rhom_tp=Rhom_tp;
options.Raux=RV;
if strcmp(options.linAlg,'wrapping-free')
    options.Rpar=interval(RV);
else
    options.Rpar=RV;
end
options.Rtrans=obj.taylor.Rtrans;

%total solution
if isa(Rinit,'mptPolytope')
    %convert zonotopes to polytopes
    Radd=mptPolytope(RV);
    Rtotal=Rhom+Radd;
    Rtotal_tp=Rhom_tp+Radd;
else
    %original computation
    Rtotal=Rhom+RV;
    Rtotal_tp=Rhom_tp+RV;
end

%write results to reachable set struct Rfirst
Rfirst.tp=Rtotal_tp;
Rfirst.ti=Rtotal;




%------------- END OF CODE --------------