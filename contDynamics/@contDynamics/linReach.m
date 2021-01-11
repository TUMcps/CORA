function [Rti,Rtp,dimForSplit,options] = linReach(obj,options,Rstart)
% linReach - computes the reachable set after linearization
%
% Syntax:  
%    [Rti,Rtp,dimForSplit,options] = linReach(obj,options,Rstart)
%
% Inputs:
%    obj - nonlinearSys or nonlinearParamSys object
%    options - struct with algorithm settings
%    Rstart - initial reachable set
%
% Outputs:
%    Rti - reachable set for time interval
%    Rtp - reachable set for time point
%    dimForSplit - dimension that is split to reduce the lin. error
%    options - struct with algorithm settings
%
% References: 
%   [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%       uncertain parameters using conservative linearization"
%   [2] M. Althoff et al. "Reachability analysis of nonlinear systems using 
%       conservative polynomialization and non-convex sets"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: initReach, post

% Author:       Matthias Althoff, Niklas Kochdumper, Mark Wetzlinger
% Written:      17-January-2008
% Last update:  29-June-2009
%               23-July-2009
%               10-July-2012
%               18-September-2012
%               09-August-2016
%               12-September-2017
%               02-January-2020 (NK, restructured the function)
%               22-April-2020 (MW, simplification)
% Last revision: ---

%------------- BEGIN CODE --------------

% extract initial set and abstraction error
Rinit = Rstart.set;
abstrerr = Rstart.error;

% linearize the nonlinear system
[obj,linSys,linOptions] = linearize(obj,options,Rinit); 

% translate Rinit by linearization point
Rdelta = Rinit + (-obj.linError.p.x);

% compute reachable set of the linearized system
if isa(obj,'nonlinParamSys') && isa(options.paramInt,'interval')
    [linSys,R] = initReach_inputDependence(linSys,Rdelta,linOptions);
else
    R = initReach(linSys,Rdelta,linOptions);

    if strcmp(options.alg,'poly')
        % precompute set of state differences
        Rdiff = deltaReach(linSys,Rdelta,linOptions);
        
        % precompute static abstraction error
        [H,Zdelta,errorStat,T,ind3,Zdelta3] = ...
            precompStatError(obj,Rdelta,options);
    end
end

% compute reachable set of the abstracted system including the
% abstraction error using the selected algorithm
if strcmp(options.alg,'linRem')
    [Rtp,Rti,perfInd] = linReach_linRem(obj,R,Rinit,Rdelta,options);
else

    % loop until the actual abstraction error is smaller than the 
    % estimated linearization error
    Rtp = R.tp; Rti = R.ti; perfIndCurr = inf; perfInd = 0;
    while perfIndCurr > 1 && perfInd <= 1
        % estimate the abstraction error 
        appliedError = 1.1*abstrerr;
        Verror = zonotope([0*appliedError,diag(appliedError)]);
        RallError = errorSolution(linSys,options,Verror); 

        % compute the abstraction error using the conservative
        % linearization approach described in [1]
        if strcmp(options.alg,'lin')

            % compute overall reachable set including linearization error
            Rmax = Rti+RallError;
            % compute linearization error
            [trueError,VerrorDyn] = abstrerr_lin(obj,options,Rmax);
            VerrorStat = [];
            
        % compute the abstraction error using the conservative
        % polynomialization approach described in [2]    
        else

            % compute overall reachable set including linearization error
            Rmax = Rdelta+zonotope(Rdiff)+RallError;
            % compute abstraction error
            [trueError,VerrorDyn,VerrorStat] = ...
                abstrerr_poly(obj,options,Rmax,Rdiff+RallError, ...
                                    H,Zdelta,errorStat,T,ind3,Zdelta3);

        end

        % compare linearization error with the maximum allowed error
        perfIndCurr = max(trueError./appliedError);    
        perfInd = max(trueError./options.maxError);

        abstrerr = trueError;
        
        % clean exit in case of set explosion
        if any(abstrerr > 1e+100)
            throw(printExplosionError());
        end
    end

    % translate reachable sets by linearization point
    Rti = Rti+obj.linError.p.x;
    Rtp = Rtp+obj.linError.p.x;

    % compute the reachable set due to the linearization error
    if isa(linSys,'linParamSys')
        Rerror = errorSolution(linSys,options,VerrorDyn);
    else
        Rerror = errorSolution(linSys,options,VerrorDyn,VerrorStat);
    end
    
    % add the abstraction error to the reachable sets
    if isa(Rerror,'polyZonotope')
        Rti=exactPlus(Rti,Rerror);
        Rtp=exactPlus(Rtp,Rerror);
    else
        Rti=Rti+Rerror;
        Rtp=Rtp+Rerror;
    end
end

% determine the best dimension to split the set in order to reduce the
% linearization error
dimForSplit = [];

if perfInd > 1
    dimForSplit = select(obj,options,Rstart);
end

% store the linearization error
Rtp_.set = Rtp;
Rtp_.error = abstrerr;
Rtp = Rtp_;


end



% Auxiliary Functions -----------------------------------------------------

function [Rtp,Rti,perfInd] = linReach_linRem(obj,R,Rinit,Rdelta,options)
% Compute the reachable set for the linearized system using an algorithm
% that is based on the linearization of the Lagrange remainder
    
    % compute the reachable set for the linearized system
    options.alg = 'lin';
    
    [obj,linSys,linOptions] = linearize(obj,options,Rinit); 
    if isa(linSys,'linParamSys')
        linOptions.compTimePoint = 1;
    end
    if isa(obj,'nonlinParamSys') && isa(options.paramInt,'interval')
        [~,Rlin] = initReach_inputDependence(linSys,Rdelta,linOptions);
    else
        Rlin = initReach(linSys,Rdelta,linOptions);
    end
    
    % compare the computed reachable set to the reachable set of the
    % linearized system in order to decide if splitting is required
    Ro_int = interval(Rlin.ti);
    Rti_int = interval(R.ti);
    
    assert(Ro_int<=Rti_int,'Bug: should be always contained');
    
    trueError = max(abs(Rti_int.inf-Ro_int.inf),abs(Rti_int.sup-Ro_int.sup));
    perfInd = max(trueError./options.maxError);
    
    % translate reachable sets by linearization point
    Rti = R.ti + obj.linError.p.x;
    Rtp = R.tp + obj.linError.p.x;
    
end

%------------- END OF CODE --------------