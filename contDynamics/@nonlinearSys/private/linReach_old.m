function [Rti,Rtp,perfInd,nr,options] = linReach(obj,options,Rinit,iter)
% linReach - computes the reachable set after linearazation and returns if
% the initial set has to be split in order to control the linearization
% error
%
% Syntax:  
%    [Rti,Rtp,split,options] = linReach(obj,options,Rinit)
%
% Inputs:
%    obj - nonlinear system object
%    options - options struct
%    Rinit - initial reachable set
%
% Outputs:
%    Rti - reachable set for time interval
%    Rti - reachable set for time point
%    split - boolean value returning if initial set has to be split
%    options - options struct to return f0
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      17-January-2008
% Last update:  29-June-2009
%               23-July-2009
%               10-July-2012
% Last revision: ---

%------------- BEGIN CODE --------------


% linearize nonlinear system
[obj,linSys,options,linOptions] = linearize(obj,options,Rinit); 

%translate Rinit by linearization point
Rdelta=Rinit+(-obj.linError.p.x);

% compute reachable set of linearized system
%ti: time interval, tp: time point
linOptions.p = obj.linError.p.x;
[linSys,R] = initReach(linSys,Rdelta,linOptions);
Rtp=R.tp;
Rti=R.ti;

perfIndCurr=inf;
while perfIndCurr>1
    
    %new technique: compute reachable set due to assumed linearization error
    %Verror = zonotope([0*obj.expFactor,diag(obj.expFactor)]);
    appliedError = 1.1*options.oldError;
    %appliedError = options.maxError;
    Verror = zonotope([0*appliedError,diag(appliedError)]);
    [RallError] = errorSolution(linSys,Verror,options); 

    %compute maximum reachable set due to maximal allowed linearization error
    Rmax=Rti+RallError;

    % obtain linearization error
    if options.advancedLinErrorComp
        %[Verror,error] = linError_quadratic(obj,options,Rmax);
        if options.tensorOrder<=2
            [Verror,error] = linError_mixed_noInt(obj,options,Rmax); 
        else
            [Verror,error] = linError_thirdOrder(obj,options,Rmax);
        end
    else
        [error] = linError(obj,options,Rmax);
        Verror=zonotope([0*error,diag(error)]);
    end

    %store old error
    options.oldError = error;

    %compute performance index of linearization error
    perfIndCurr = max(error./appliedError)
    
    if perfIndCurr > 1
        disp('investigate');
    end
end
% compute reachable set due to the linearization error
if strcmp('quadZonotope',class(Verror))
    [Rerror] = errorSolutionQuad(linSys,Verror,options); %<-- this has to change!
else
    [Rerror] = errorSolution(linSys,Verror,options); %<-- this has to change!
end

%translate reachable sets by linearization point
Rti=Rti+obj.linError.p.x;
Rtp=Rtp+obj.linError.p.x;

perfInd = max(error./options.maxError)

if (perfInd>1) && (iter==1)
    %find best split
    nr=select(obj,options,Rinit,obj.expFactor);
else
    nr=[];
end

%add intervalhull of actual error
if strcmp('quadZonotope',class(Rerror))
    Rti=exactPlus(Rti,Rerror);
    Rtp=exactPlus(Rtp,Rerror);
else
    Rti=Rti+Rerror;
    Rtp=Rtp+Rerror;
end


%------------- END OF CODE --------------