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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       15-July-2019 (from @linearSys > post_Euclidean.m)
% Last update:   16-February-2021 (MW, update 'fromStart')
%                16-November-2021 (MW, include disturbance set W)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load data

eAt = obj.taylor.eAt;
eAtInt = obj.taylor.eAtInt;
F = obj.taylor.F;
inputF = obj.taylor.inputF;

if isfield(options,'W')
    Wcenter = center(options.W);
else
    Wcenter = 0;
end

% solution due to constant inputs
vTrans = obj.B*options.uTrans + Wcenter;

if ~isempty(obj.c)
    vTrans = vTrans + obj.c;
end

Rtrans = eAtInt*zonotope(vTrans);
% effect should only be considered once in a single time interval:
inputCorr = inputF*zonotope(vTrans);

% homogeneous time-point solution
Rinit = options.Rhom_tp;
if strcmp(options.linAlg,'wrapping-free')
    Rhom_tp = eAt*Rinit + center(Rtrans);
elseif any(strcmp(options.linAlg,{'standard','fromStart'}))
    Rhom_tp = eAt*Rinit + Rtrans;
end

% homogeneous time-interval solution
% note: inputCorr considered later
if isa(Rinit,'polyZonotope') 
    Rhom = enclose(Rinit,Rhom_tp) + F*zonotope(Rinit);
elseif isa(Rinit,'zonoBundle') 
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit.Z{1};
else
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit;
end


end

% ------------------------------ END OF CODE ------------------------------
