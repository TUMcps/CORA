function [Rti,Rtp,Rti_y,perfInd,nr,options] = linReach(obj,options,Rinit,Rinit_y,iter)
% linReach - computes the reachable set after linearazation and returns if
% the initial set has to be split in order to control the linearization
% error
%
% Syntax:  
%    [Rti,Rtp,perfInd,nr,options] = linReach(obj,options,Rinit,iter)
%
% Inputs:
%    obj - nonlinear DAE system object
%    options - options struct
%    Rinit - initial reachable set
%    iter - flag for activating iteration
%
% Outputs:
%    Rti - reachable set for time interval
%    Rti - reachable set for time point
%    perfInd - performance index
%    nr - number of generator that should be split
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
% Written:      21-November-2011
% Last update:  28-May-2013
% Last revision:---

%------------- BEGIN CODE --------------

% linearize nonlinear system
[obj,linSys,options,linOptions] = linearize(obj,options,Rinit,Rinit_y); 

%translate Rinit by linearization point
Rdelta = Rinit + (-obj.linError.p.x);

% compute reachable set of linearized system
%ti: time interval, tp: time point
linOptions.p = obj.linError.p.x;
R = initReach(linSys,Rdelta,linOptions);
Rtp=R.tp;
Rti=R.ti;

perfIndCurr_x=inf;
perfIndCurr_y=inf;
perfInd=0;
expFactor = 1.1;
while ((perfIndCurr_x>1) || (perfIndCurr_y>1)) && (perfInd<=1) 

    %new technique: compute reachable set due to assumed linearization error
    %Verror = zonotope([0*obj.expFactor,diag(obj.expFactor)]);
    if perfIndCurr_x>1
        appliedError_x = expFactor*options.oldError_x;
    end
    if perfIndCurr_y>1
        appliedError_y = expFactor*options.oldError_y;
    end

    %convert error to zonotope
    Verror_x = zonotope([0*appliedError_x,diag(appliedError_x)]);
    Verror_y = zonotope([0*appliedError_y,diag(appliedError_y)]);
    Verror = Verror_x + obj.linError.CF_inv*Verror_y;

    [RallError] = errorSolution(linSys,Verror,options); 

    %compute maximum reachable set due to maximal allowed linearization error
    Rmax = Rti + RallError;

    % obtain linearization error
    if options.advancedLinErrorComp == 1
        [Verror, error, error_x, error_y, Rti_y] = linError_mixed_noInt(obj, options, Rmax, Verror_y); 
    elseif options.advancedLinErrorComp == 2
        [Verror, error, error_x, error_y, Rti_y] = linError_thirdOrder(obj, options, Rmax, Verror_y);
    elseif options.advancedLinErrorComp == 21
%         if strcmp(options.category,'powerSystem')
%             %[Verror, error, error_x, error_y, Rti_y] = linError_mixed_noInt_comp(obj, options, Rmax, Verror_y);
%             [Verror, error, error_x, error_y, Rti_y] = linError_powDyn_mixed_noInt_comp(obj, options, Rmax, Verror_y);
%         else
            [Verror, error, error_x, error_y, Rti_y] = linError_mixed_noInt_comp(obj, options, Rmax, Verror_y);
            %[Verror, error, error_x, error_y, Rti_y] = linError_mixed_noInt_comp_parallel(obj, options, Rmax, Verror_y);
%         end
    elseif options.advancedLinErrorComp == 22
        [Verror, error, error_x, error_y, Rti_y] = linError_thirdOrder_comp(obj, options, Rmax, Verror_y);
    elseif options.advancedLinErrorComp == 31
        [Verror, error, error_x, error_y, Rti_y] = linError_mixed_specialTOP(obj, options, Rmax, Verror_y); 
    else
        [Verror, error, error_x, error_y, Rti_y] = linError(obj, options, Rmax, Verror_y);
    end

    %store old error
    options.oldError = error;
    options.oldError_x = error_x;
    options.oldError_y = error_y;

    %compute performance index of linearization error
    perfIndCurr_x = max(error_x./appliedError_x)
    perfInd_x = max(error_x./options.maxError_x)
    perfIndCurr_y = max(error_y./appliedError_y)
    perfInd_y = max(error_y./options.maxError_y)
    
    if (perfIndCurr_x > 1) || (perfIndCurr_y > 1)
        disp('investigate');
    end

    perfInd = max(perfInd_x, perfInd_y);
end

% compute reachable set due to the linearization error
[Rerror] = errorSolution(linSys,Verror,options);

%translate reachable sets by linearization point
Rti=Rti+obj.linError.p.x;
Rtp=Rtp+obj.linError.p.x;

if perfInd>0.8
    disp('investigate');
end

if (perfInd>1) && (iter==1)
    %find best split
clc    nr=select(obj,options,Rinit,obj.expFactor);
else
    nr=[];
end

%add interval of actual error
Rti=Rti+Rerror;
Rtp=Rtp+Rerror;


%------------- END OF CODE --------------