function [R,res] = reach(obj,params,varargin)
% reach - computes the reachable set for linear systems
%
% Syntax:  
%    R = reach(obj,params)
%    R = reach(obj,params,options)
%    [R,res] = reach(obj,param,spec)
%    [R,res] = reach(obj,param,options,spec)
%
% Inputs:
%    obj - continuous system object
%    params - model parameters
%    options - options for the computation of reachable sets
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res  - 1 if specifications are satisfied, 0 if not
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Mark Wetzlinger
% Written:       26-June-2019
% Last update:   08-Oct-2019
%                23-April-2020 (added params)
% Last revision: ---


%------------- BEGIN CODE --------------

    % parse input arguments
    spec = [];
    if nargin == 2
        options.linAlg = 'adap'; 
    elseif nargin == 3
        if isa(varargin{1},'specification')
            spec = varargin{1};
            options.linAlg = 'adap';
        else
            options = varargin{1};
        end
    else
        options = varargin{1};
        spec = varargin{2};
    end
    
    % options preprocessing
    options = params2options(params, options);
    options = checkOptionsReach(obj,options);
    
    if ~isempty(spec)
    	options.specification = spec; 
    end

    % decide which reach function to execute by options.linAlg
    if strcmp(options.linAlg,'adap')
        [Rout,Rout_tp,res,deltaT] = reach_adaptive(obj, options);
        tVec = cumsum(deltaT) + options.tStart;
    else
        % all below, const. time step sizes
        tVec = (options.tStart+options.timeStep:options.timeStep:options.tFinal)';
        if strcmp(options.linAlg,'standard')
            [Rout,Rout_tp,res] = reach_standard(obj, options);
        elseif strcmp(options.linAlg,'wrapping-free')
            [Rout,Rout_tp,res] = reach_wrappingfree(obj, options);
        elseif strcmp(options.linAlg,'fromStart')
            [Rout,Rout_tp,res] = reach_fromStart(obj, options);
        elseif strcmp(options.linAlg,'decomp')
            [Rout,Rout_tp,res] = reach_decomp(obj, options);
        elseif strcmp(options.linAlg,'krylov')
            [Rout,Rout_tp,res] = reach_krylov(obj, options);
        end
    end

    % create object of class reachSet
    R = createReachSetObject(Rout,Rout_tp,tVec,options);

end


% Auxiliary Functions -----------------------------------------------------

function R = createReachSetObject(Rout,Rout_tp,tVec,options)
% create and object of class reachSet that stores the reachable set

    % time-point reachable set
    timePoint.set = Rout_tp;
    timePoint.time = num2cell(tVec);
    
    if length(timePoint.time) ~= length(timePoint.set)
       timePoint.time = timePoint.time(1:length(timePoint.set)); 
    end
    
    % time interval reachable set
    timeInt.set = Rout;
    timeInt.time = cell(length(Rout),1);
    
    t = options.tStart;
    
    for i = 1:length(Rout)
        timeInt.time{i} = interval(t,tVec(i));
        t = tVec(i);
    end
    
    % construct object of class reachSet
    if strcmp(options.linAlg,'decomp')
        R = [];
        for i = 1:length(timePoint.set{1})
           timePoint_ = timePoint;
           timePoint_.set = cellfun(@(x) x{i},timePoint.set,'UniformOutput',false);
           timeInt_ = timeInt;
           timeInt_.set = cellfun(@(x) x{i},timeInt.set,'UniformOutput',false);
           if ~isempty(timeInt_.set{1})
               if isempty(R)
                   R = reachSet(timePoint_,timeInt_); 
               else
                   R = add(R,reachSet(timePoint_,timeInt_));
               end
           end
        end
    else
        R = reachSet(timePoint,timeInt);
    end
end

%------------- END OF CODE --------------