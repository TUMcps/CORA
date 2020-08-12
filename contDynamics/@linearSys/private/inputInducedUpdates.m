function [Rhom,Rhom_tp,Rtrans,inputCorr] = inputInducedUpdates(obj,options)
% inputInducedUpdates - recalculates set if options.uTrans has changed
%
% Syntax:  
%    [Rhom,Rhom_tp,Rtrans,inputCorr] = inputInducedUpdates(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    Rhom - homogeneous time interval solution
%    Rhom_tp - homogeneous time point solution
%    Rtrans - particular solution using uTrans
%    inputCorr - bloating due to input change (only time interval)
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       15-July-2019 (from @linearSys > post_Euclidean.m)
% Last update:   ---
% Last revision: ---


%------------- BEGIN CODE --------------

% load data
eAt = obj.taylor.eAt;
eAtInt = obj.taylor.eAtInt;
F = obj.taylor.F;
inputF = obj.taylor.inputF;

% solution due to constant inputs
vTrans = obj.B*options.uTrans;
Rtrans = eAtInt*zonotope(vTrans);
% effect should only be considered once in a single time interval:
inputCorr = inputF*zonotope(vTrans);

% homogeneous time point solution
if strcmp(options.linAlg,'fromStart')
    Rinit = options.R0;
    Rhom_tp = eAt*Rinit + Rtrans;
else
    Rinit = options.Rhom_tp;
    if strcmp(options.linAlg,'wrapping-free')
        Rhom_tp = eAt*Rinit + center(Rtrans);
    elseif strcmp(options.linAlg,'standard')
        Rhom_tp = eAt*Rinit + Rtrans;
    end
end


% homogeneous time interval solution
% note: inputCorr considered later
if isa(Rinit,'quadZonotope') 
    Rhom = enclose(Rinit,Rhom_tp) + F*zonotope(Rinit);
elseif isa(Rinit,'zonoBundle') 
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit.Z{1};
else
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit;
end


end

%------------- END OF CODE --------------