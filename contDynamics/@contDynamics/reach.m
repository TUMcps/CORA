function [R,res] = reach(obj,params,options,varargin)
% reach - computes the reachable continuous set for the entire time horizon
%         of a continuous system
%
% Syntax:  
%    R = reach(obj,params,options)
%    [R,res] = reach(obj,params,options,spec)
%
% Inputs:
%    obj - continuous system object
%    params - parameter defining the reachability problem
%    options - options for the computation of reachable sets
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the computed reachable set
%    res  - 1 if specifications are satisfied, 0 if not
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      08-August-2016
% Last update:  22-September-2016
%               28-July-2017 By Elguindy - the reachable set of alg.
%                            variables lines 64-68
%                            Additional output Rcont_y is now included
%               20-March-2018 (NK, output sets as additional output)
%               19-May-2020 (MW, error handling for exploding sets)
% Last revision:---

%------------- BEGIN CODE --------------

    res = 1;

    % options preprocessing
    options = params2options(params,options);
    options = checkOptionsReach(obj,options);
    
    spec = [];
    if nargin >= 4
       spec = varargin{1}; 
    end

    % compute symbolic derivatives
    if isa(obj,'nonlinearSys') || isa(obj,'nonlinDASys') || isa(obj,'nonlinParamSys')
       derivatives(obj,options); 
    end

    % obtain factors for initial state and input solution time step
    r = options.timeStep;
    for i = 1:(options.taylorTerms+1)  
        options.factor(i) = r^(i)/factorial(i);    
    end

    % if a trajectory should be tracked
    if isfield(options,'uTransVec')
        options.uTrans = options.uTransVec(:,1);
    end
    options.t = options.tStart;
    
    % time period
    tVec = options.tStart:options.timeStep:options.tFinal;

    % initialize cell-arrays that store the reachable set
    timeInt.set = cell(length(tVec)-1,1);
    timeInt.time = cell(length(tVec)-1,1);
    timePoint.set = cell(length(tVec)-1,1);
    timePoint.time = cell(length(tVec)-1,1);
    if isa(obj,'nonlinDASys')
        timeInt.algebraic = cell(length(tVec)-1,1);
    end

    % initialize reachable set computations
    try
        [Rnext, options] = initReach(obj, options.R0, options);
    catch ME
        % if error from set explosion, return corresponding information
        R = [];
        reportReachError(ME,options.tStart,1);
        return
    end
    
    % loop over all reachability steps
    for i = 2:length(tVec)-1
        
        % save reachable set in cell structure
        timeInt.set{i-1} = Rnext.ti; 
        timeInt.time{i-1} = interval(tVec(i-1),tVec(i));
        timePoint.set{i-1} = Rnext.tp;
        timePoint.time{i-1} = tVec(i);
        
        if isa(obj,'nonlinDASys')
            timeInt.algebraic{i-1}  = Rnext.y;
        end
        
        % check specification
        if ~isempty(spec)
           if ~check(spec,Rnext.ti)
               res = 0;
               R = createReachSetObject(timeInt,timePoint);
               return;
           end
        end

        % increment time and set counter
        t = tVec(i);
        options.t = t;
        if isfield(options,'verbose') && options.verbose 
            disp(t);
        end

        % if a trajectory should be tracked
        if isfield(options,'uTransVec')
            options.uTrans = options.uTransVec(:,i);
        end

        % compute next reachable set
        try
            [Rnext,options] = post(obj,Rnext,options);
        catch ME
            % if error from set explosion, return corresponding information
            R = createReachSetObject(timeInt,timePoint);
            reportReachError(ME,t,i);
            return
        end
        
    end
    
    % check specification
    if ~isempty(spec)
       if ~check(spec,Rnext.ti)
           res = 0;
       end
    end

    % save last reachable set in cell structure
    timeInt.set{end} = Rnext.ti; 
    timeInt.time{end} = interval(tVec(end-1),tVec(end));
    timePoint.set{end} = Rnext.tp; 
    timePoint.time{end} = tVec(end);

    if isfield(Rnext,'y')
        timeInt.algebraic{end} = Rnext.y;
    end
    
    % construct reachset object
    R = createReachSetObject(timeInt,timePoint);
end


% Auxiliary Functions -----------------------------------------------------

function R = createReachSetObject(timeInt,timePoint)
% create and object of class reachSet that stores the reachable set

    % remove empty cells
    if isempty(timeInt.set{end})
       ind = ~cellfun('isempty',timeInt.set);
       timePoint.set = timePoint.set(ind);
       timePoint.time = timePoint.time(ind);
       timeInt.set = timeInt.set(ind);
       timeInt.time = timeInt.time(ind);
       if isfield(timeInt,'algebraic')
           timeInt.algebraic = timeInt.algebraic(ind);
       end
    end

    % check if sets are split is required
    if ~iscell(timeInt.set{1})
        if ~isempty(timePoint.set{1}) || isa(timePoint.set{1},'zonoBundle')
            R = reachSet(timePoint,timeInt);
        else
            R = reachSet([],timeInt); 
        end
    else

        % no splitting occured -> copy sets
        if length(timeInt.set{end}) == 1   
            timeInt.set = cellfun(@(x) x{1},timeInt.set,'UniformOutput',false);
            if isfield(timeInt,'algebraic')
                timePoint.set = cellfun(@(x) x{1}.set,timePoint.set,'UniformOutput',false);
                timeInt.algebraic = cellfun(@(x) x{1},timeInt.algebraic,'UniformOutput',false);
            else
                timePoint.set = cellfun(@(x) x{1}.set,timePoint.set,'UniformOutput',false);
            end
            R = reachSet(timePoint,timeInt);

        % splitted sets -> bring to correct format
        else

            ind = 1:length(timePoint.set{1});
            timeInt_ = {}; timePoint_ = {};
            parent = zeros(length(timePoint.set{1}),1);

            % loop over all time steps
            for i = 1:length(timePoint.set)

               % modify index vector
               ind_ = zeros(length(timePoint.set{i}),1);
               maxInd = max(ind);

               for j = 1:length(timePoint.set{i})
                  if isfield(timePoint.set{i}{j},'parent') && i > 1
                      ind_(j) = maxInd + 1;
                      parent = [parent; ind(timePoint.set{i}{j}.parent)];
                      maxInd = maxInd + 1;
                  else
                      ind_(j) = ind(timePoint.set{i}{j}.prev);
                  end
               end
               ind = ind_;

               % copy entries
               for j = 1:length(timePoint.set{i})
                   if ind(j) <= length(timeInt_)
                      timeInt_{ind(j)}.set = [timeInt_{ind(j)}.set; timeInt.set{i}(j)];
                      timeInt_{ind(j)}.time = [timeInt_{ind(j)}.time; timeInt.time(i)];
                      if isfield(timeInt,'algebraic')
                          timeInt_{ind(j)}.algebraic = [timeInt_{ind(j)}.algebraic; timeInt.algebraic{i}(j)];
                      end
                      timePoint_{ind(j)}.set = [timePoint_{ind(j)}.set; {timePoint.set{i}{j}.set}];
                      timePoint_{ind(j)}.time = [timePoint_{ind(j)}.time; timePoint.time(i)];
                   else
                      timeInt_{ind(j)}.set = timeInt.set{i}(j);
                      timeInt_{ind(j)}.time = timeInt.time(i);
                      if isfield(timeInt,'algebraic')
                          timeInt_{ind(j)}.algebraic = timeInt.algebraic{i}(j);
                      end
                      timePoint_{ind(j)}.set = {timePoint.set{i}{j}.set};
                      timePoint_{ind(j)}.time = timePoint.time(i);
                   end
               end
            end

            % generate reachSet object
            for i = 1:length(timeInt_)
                
                R_ = reachSet(timePoint_{i},timeInt_{i});

                if i == 1
                    R = R_; 
                else
                    if parent(i) > 0
                        R = add(R,R_,parent(i));
                    else
                        R = add(R,R_);
                    end
                end
            end
        end
    end
end

%------------- END OF CODE --------------