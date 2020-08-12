function [Rnext,options] = post_Euclidean(obj,options)
% post_Euclidean - computes the reachable continuous set for one time step
% in the untransformed space
%
% Syntax:  
%    [Rnext,options] = post_Euclidean(obj,options)
%
% Inputs:
%    obj - linearSys object
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rnext - reachable set of the next time step
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      07-May-2007 
% Last update:  27-April-2009
%               29-June-2009
%               08-August-2011
%               25-July-2016 (intervalhull replaced by interval)
%               07-July-2018 (wrapping-free approach is no longer the default option)
%               27-July-2018 (changing input updated)
%               18-September-2018 (changing input again updated)
%               07-November-2018 (arguments of function updated)
% Last revision:---

%------------- BEGIN CODE --------------

%load data
eAt=obj.taylor.eAt;

% check whether standard or wrapping-free approach should be used
if strcmp(options.linAlg,'wrapping-free')
    % option 1 (wrapping free)
    % method implemented from Algorithm 2 in
    % A. Girard, C. Le Guernic, and O. Maler, “Efficient computation of 
    % reachable sets of linear time-invariant systems with inputs,�? in 
    % Hybrid Systems: Computation and Control, ser. LNCS 3927. Springer, 
    % 2006, pp. 257–271.
    if isfield(options,'uTransVec') % check whether input vector changes
        [Rhom,Rhom_tp,Rtrans,inputCorr] = inputInducedUpdates(obj,options);
    else
        Rtrans = options.Rtrans;
        inputCorr = 0;
        Rhom=eAt*options.Rhom + center(Rtrans);
        if options.compTimePoint
            Rhom_tp=eAt*options.Rhom_tp + center(Rtrans);
        end
    end
    Raux=eAt*options.Raux;
    Rpar=options.Rpar + interval(Raux) + interval(Rtrans) + (-center(Rtrans));
elseif strcmp(options.linAlg,'standard')
    % option 2 (not wrapping-free)
    % method implemented from Algorithm 1 in
    % A. Girard, “Reachability of uncertain linear systems using 
    % zonotopes,�? in Hybrid Systems: Computation and Control, 
    % ser. LNCS 3414. Springer, 2005, pp. 291–305.
    if isfield(options,'uTransVec') % check whether input vector changes
        [Rhom,Rhom_tp,~,inputCorr] = inputInducedUpdates(obj,options);
    else
        Rtrans = options.Rtrans;
        inputCorr = 0;
        Rhom=eAt*options.Rhom + Rtrans;
        if options.compTimePoint
            Rhom_tp=eAt*options.Rhom_tp + Rtrans;
        end
    end
    Raux=eAt*options.Raux;
    Rpar=reduce(options.Rpar + Raux,options.reductionTechnique,options.zonotopeOrder);
    Rhom=reduce(Rhom,options.reductionTechnique,options.zonotopeOrder);
    if options.compTimePoint
        Rhom_tp=reduce(Rhom_tp,options.reductionTechnique,options.zonotopeOrder);
    end
else
    error("Don't use this function with the given options.linAlg");
end

%save homogeneous and particulate solution to options struct
options.Rhom=Rhom;
if options.compTimePoint
    options.Rhom_tp=Rhom_tp;
end
options.Rpar=Rpar;
options.Raux=Raux;

%write results to reachable set struct Rnext
if isa(Rhom,'mptPolytope')
    Rnext.ti=Rhom+mptPolytope(Rpar)+mptPolytope(inputCorr);
    if options.compTimePoint
        Rnext.tp=Rhom_tp+mptPolytope(Rpar);
    else
        Rnext.tp = [];
    end
else
    Rnext.ti=Rhom+zonotope(Rpar)+inputCorr;
    if options.compTimePoint
        Rnext.tp=Rhom_tp+zonotope(Rpar);
    else
        Rnext.tp = [];
    end
end

end



%------------- END OF CODE --------------