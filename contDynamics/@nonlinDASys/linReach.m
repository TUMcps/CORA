function [Rti,Rtp,Rti_y,perfInd,dimForSplit,options] = linReach(nlnsysDA,Rinit,Rinit_y,params,options,iter)
% linReach - computes the reachable set after linearazation and returns if
%    the initial set has to be split in order to control the linearization
%    error
%
% Syntax:
%    [Rti,Rtp,perfInd,dimForSplit,options] = ...
%       linReach(nlnsysDA,Rinit,Rinit_y,params,options,iter)
%
% Inputs:
%    nlnsysDA - nonlinDASys object
%    Rinit - initial reachable set (diff. variables)
%    Rinit_y - initial reachable set (alg. variables)
%    params - model parameters
%    options - options struct
%    iter - flag for activating iteration
%
% Outputs:
%    Rti - reachable set for time interval (diff. variables)
%    Rtp - reachable set for time point (diff. variables)
%    Rti_y - reachable set for time interval (alg. variables)
%    perfInd - performance index
%    dimForSplit - number of generator that should be split
%    options - options struct to return f0
%
% Example:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       21-November-2011
% Last update:   28-May-2013
%                19-May-2020 (MW, error handling for exploding sets)
%                11-January-2021 (MW, syntax change for reachSet init)
%                26-May-2022 (MA, explicit selection of linearization error computation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% extract start set and abstraction error
abstrerr_x = Rinit.error_x;
abstrerr_y = Rinit.error_y;
Rinit = Rinit.set;

% linearize nonlinear system
[nlnsysDA,linsys,linParams,linOptions] = priv_linearize(nlnsysDA,Rinit,Rinit_y,params,options); 

%translate Rinit by linearization point
Rdelta = Rinit + (-nlnsysDA.linError.p.x);

% compute reachable set of linearized system
% ti: time interval, tp: time point
linOptions.p = nlnsysDA.linError.p.x;
[Rtp,Rti] = oneStep(linsys,Rdelta,linParams.U,linParams.uTrans,options.timeStep,options.taylorTerms);

% performance indices for splitting
perfIndCurr_x = Inf; perfIndCurr_y = Inf;
perfInd = 0;
expFactor = 1.1;

% loop until assumed linearization error larger than computed a posteriori
while ((perfIndCurr_x > 1) || (perfIndCurr_y > 1)) && (perfInd <= 1) 

    if perfIndCurr_x > 1
        appliedError_x = expFactor * abstrerr_x;
    end
    if perfIndCurr_y > 1
        appliedError_y = expFactor * abstrerr_y;
    end

    %convert error to zonotope
    Verror_x = zonotope(0*appliedError_x,diag(appliedError_x));
    Verror_y = zonotope(0*appliedError_y,diag(appliedError_y));
    Verror = Verror_x + nlnsysDA.linError.CF_inv * Verror_y;

    RallError = particularSolution_timeVarying(linsys,Verror,...
        options.timeStep,options.taylorTerms); 

    %compute maximum reachable set due to maximal allowed linearization error
    Rmax = Rti + RallError;

    % obtain linearization error
    if options.tensorOrder == 2
%         [Verror, error, error_x, error_y, Rti_y] = ...
%            priv_linError(nlnsysDA, options, Rmax, Verror_y);
        if ~isfield(options,'index')
            % conventional computation
            [Verror, error, error_x, error_y, Rti_y] = ...
                priv_linError_mixed_noInt(nlnsysDA, Rmax, Verror_y, params, options);
        else
            % compositional computation
            [Verror, error, error_x, error_y, Rti_y] = ...
                priv_linError_mixed_noInt_comp(nlnsysDA, Rmax, Verror_y, params, options);
        end
    elseif options.tensorOrder == 3
        [Verror, error, error_x, error_y, Rti_y] = ...
            priv_linError_thirdOrder(nlnsysDA, Rmax, Verror_y, params, options);
    end
    
    
    %compute performance index of linearization error
    perfIndCurr_x = max(error_x ./ appliedError_x);
    perfInd_x = max(error_x ./ options.maxError_x);
    perfIndCurr_y = max(error_y ./ appliedError_y);
    perfInd_y = max(error_y ./ options.maxError_y);

    % compute overall performance index
    perfInd = max(perfInd_x, perfInd_y);
    
    % store error
    abstrerr_x = error_x;
    abstrerr_y = error_y;
    
    % clean exit in case of set explosion
    if any(error > 1e+100)
        throw(CORAerror('CORA:reachSetExplosion'));
    end
end

% compute reachable set due to the linearization error
Rerror = particularSolution_timeVarying(linsys,Verror,...
    options.timeStep,options.taylorTerms);

%translate reachable sets by linearization point
Rti = Rti + nlnsysDA.linError.p.x;
Rtp = Rtp + nlnsysDA.linError.p.x;

if perfInd > 0.8
    disp('investigate');
end

dimForSplit = [];
if perfInd > 1 && iter == 1
    % find best split
    dimForSplit = priv_select(nlnsysDA,Rinit,Rinit_y,params,options,iter);
end

%add interval of actual error
Rti = Rti + Rerror;
Rtp = Rtp + Rerror;

% store the linearization error, correct syntax for reachSet object later
Rtp_.set = Rtp;
Rtp_.error_x = error_x;
Rtp_.error_y = error_y;
Rtp = Rtp_;

% ------------------------------ END OF CODE ------------------------------
