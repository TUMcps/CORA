function [Rti,RtpEnh,nr,perfInd] = linReach(obj,options,Rstart,recur)
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
%    recur - flag if recursive calls are allowed
%
% Outputs:
%    Rti - reachable set for time interval
%    Rti - reachable set for time point
%    split - boolean value returning if initial set has to be split
%    perfInd - performance index

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
%               18-September-2012
%               09-August-2016
%               12-September-2017
% Last revision: ---

%------------- BEGIN CODE --------------

%init nr for splitting
nr = [];

%extract set and error
Rinit = Rstart.set;
error = Rstart.error;

% linearize nonlinear system
[obj,linSys,linOptions] = linearize(obj,options,Rinit); 

%translate Rinit by linearization point
Rdelta=Rinit+(-obj.linError.p.x);

% compute reachable set of linearized system
%ti: time interval, tp: time point
linOptions.p = obj.linError.p.x;
R = initReach(linSys,Rdelta,linOptions);
if options.advancedLinErrorComp
    Rdiff = deltaReach(linSys,Rdelta,linOptions);
end

Rtp=R.tp;
Rti=R.ti;

if options.advancedLinErrorComp == 2 % no error considered
    Rerror = zonotope(0);
else
    perfIndCurr=inf;
    perfInd=0;
    while (perfIndCurr>1) && (perfInd<=1) 

        %new technique: compute reachable set due to assumed linearization error
        appliedError = 1.1*error;
        Verror = zonotope([0*appliedError,diag(appliedError)]);
        [RallError] = errorSolution(linSys,Verror,options); 

        % obtain linearization error
        if options.advancedLinErrorComp
            %[Verror,error] = linError_quadratic(obj,options,Rmax);
            
            
            %compute maximum reachable set due to maximal allowed linearization error
            %Instead of Rmax=Rtp+Rdiff+RallError, one can also compute
            %Rmax=Rti+RallError; the latter approach is more accurate,
            %but results in a zonotope with more generators. Since the
            %number of generators has to be reduced in any case for the
            %error computation, the solution which already has less
            %generators is preferred.
            % Althof: "Reachability Analysis of Nonlinear Systems using
            % Conservative Polynomialization and Non-Convex Sets" 
            % Algorithm 1, line 7
            try
                Rmax=Rdelta+zonotope(Rdiff)+RallError;
            catch
                Rmax=Rdelta+Rdiff+RallError; 
            end
 
            % compute linearization error
            if options.tensorOrder<=2

                [VerrorDyn,trueError] = linError_mixed_noInt(obj,options,Rmax); 
                VerrorStat = [];
            else
                [VerrorStat,VerrorDyn,trueError] = linError_thirdOrder(obj,options,Rmax,Rdelta,Rdiff+RallError);
            end
        else
            %compute maximum reachable set due to maximal allowed linearization error
            Rmax=Rti+RallError;         

            if options.tensorOrder <= 2
                [trueError] = linError(obj,options,Rmax);
                VerrorDyn=zonotope([0*trueError,diag(trueError)]);
                VerrorStat = [];
            else
                [VerrorDyn,trueError] = linError_higherOrder(obj,Rmax,options);
                VerrorStat = [];
            end
        end

        %store old error
        error = trueError;

        %compute performance index of linearization error
        perfIndCurr = max(trueError./appliedError);    
        perfInd = max(trueError./options.maxError);
    end
    % compute reachable set due to the linearization error
    if ~isempty(VerrorStat)
        [Rerror] = errorSolutionQuad(linSys,VerrorStat,VerrorDyn,options); 
    else
        [Rerror] = errorSolution(linSys,VerrorDyn,options);
    end
    
    if (perfInd>1)  && recur
        %find best split
        nr=select(obj,options,Rstart);
    end
end

%translate reachable sets by linearization point
Rti=Rti+obj.linError.p.x;
Rtp=Rtp+obj.linError.p.x;


%add intervalhull of actual error
if isa(Rerror,'quadZonotope')
    Rti=exactPlus(Rti,Rerror);
    Rtp=exactPlus(Rtp,Rerror);
else
    Rti=Rti+Rerror;
    Rtp=Rtp+Rerror;
end

%add error to set
RtpEnh.set=Rtp;
RtpEnh.error=error;

%------------- END OF CODE --------------