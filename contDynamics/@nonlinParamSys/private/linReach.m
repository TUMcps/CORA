function [Rti,RtpEnh,nr,options] = linReach(obj,options,Rstart,recur)
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

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      17-January-2008
% Last update:  29-June-2009
%               23-July-2009
%               16-August-2016
%               07-June-2017
%               16-July-2017 (NK)
% Last revision: ---

%------------- BEGIN CODE --------------

%extract set and error
Rinit = Rstart.set;
Error = Rstart.error;

% linearize nonlinear system
[obj,linSys,linOptions] = linearize(obj,options,Rinit); 

%translate Rinit by linearization point
Rdelta=Rinit+(-obj.linError.p.x);

% compute reachable set of linearized system
%ti: time interval, tp: time point
if isa(options.paramInt,'interval')
    [linSys,R] = initReach_inputDependence(linSys,Rdelta,linOptions);
elseif isnumeric(options.paramInt)
    R = initReach(linSys,Rdelta,linOptions);
    if options.advancedLinErrorComp
        Rdiff = deltaReach(linSys,Rdelta,linOptions);
    end
end
Rtp=R.tp;
Rti=R.ti;

perfIndCurr=inf;
perfInd=0;
while (perfIndCurr>1) && (perfInd<=1) 
    
    %new technique: compute reachable set due to assumed linearization error
    appliedError = 1.1*Error;
    Verror = zonotope([0*appliedError,diag(appliedError)]);
    [RallError] = errorSolution(linSys,Verror,options);     

    % obtain linearization error
    if options.advancedLinErrorComp

        %ONLY IMPLEMENTED FOR CONSTANT PARAMETERS:
        if isa(options.paramInt,'interval')
           error('Not implemented. So far "option.advancedLinErrorComp" is only implemented for constant parameters.'); 
        end      
        
        % compute maximum reachable set due to maximal allowed linearization error
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
            %ONLY IMPLEMENTED FOR CONSTANT PARAMETERS:
            if isa(options.paramInt,'interval')
               error('Not implemented. So far options.tensorOrder > 2 is only implemented for constant parameters.'); 
            end   
            
            [VerrorDyn,trueError] = linError_higherOrder(obj,Rmax,options);
            VerrorStat = [];
        end
    end

    %store old error
    Error = trueError;

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

%translate reachable sets by linearization point
Rti=Rti+obj.linError.p.x;
Rtp=Rtp+obj.linError.p.x;

if (perfInd>1)  && recur
    %find best split
    nr=select(obj,options,Rstart);
elseif recur==0
    %return performance index instead of number if function is called
    %from select() using recur=0
    nr=perfInd;
else
    nr=[];
end

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
RtpEnh.error=Error;

%------------- END OF CODE --------------