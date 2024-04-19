function res = find(simRes,prop,val)
% find - get simRes object that satisfy a certain condition
%
% Syntax:
%    res = find(simRes,prop,val)
%
% Inputs:
%    simRes - simResult object
%    prop - property for condition ('location', 'time')
%    val - value for property
%
% Outputs:
%    res - all parts of simRes trajectories that satisfy the condition
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet/find

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       16-May-2023
% Last update:   10-April-2024 (TL, added 'time')
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments (only first two)
inputArgsCheck({{simRes,'att','simResult'};
                {prop,'str',{'location', 'time'}}});

% check which property desired
switch prop

    case 'location'
        % get all simRes trajectories with a given location number

        % quick check: all trajectories from continuous-time system?
        if isscalar(val) && val == 0 ...
                && all(arrayfun(@(x) all(x.loc==val),simRes,'UniformOutput',true))
            res = simRes; return
        end

        res = [];

        % loop over all trajectories
        for r=1:length(simRes)

            % init properties of simResult object
            x = {}; y = {}; a = {}; t = {}; loc = {};

            % loop over all individual parts
            for part=1:length(simRes(r).t)

                % check which parts satisfy loc property
                if all(simRes(r).loc(part,:) == val')
                    % save location
                    loc{end+1,1} = val';
    
                    % save corresponding time
                    t{end+1,1} = simRes(r).t{part};
    
                    % save state vector (if given)
                    if ~isempty(simRes(r).x)
                        x{end+1,1} = simRes(r).x{part};
                    end
                    % save output vector (if given)
                    if ~isempty(simRes(r).y)
                        y{end+1,1} = simRes(r).y{part};
                    end
                    % save algebraic vector (if given)
                    if ~isempty(simRes(r).a)
                        a{end+1,1} = simRes(r).a{part};
                    end
                end

            end

            % append resulting partial simRes object (only if non-empty)
            if ~isempty(t)
                res = [res; simResult(x,t,loc,y,a)];
            end
        end

    case 'time'

        % read val
        if ~isa(val,'interval')
            val = interval(val-1e-10,val+1e-10);
        end

        res = [];

        % loop over all trajectories
        for r=1:length(simRes)

            % init properties of simResult object
            x = {}; y = {}; a = {}; t = {}; loc = {};

            % loop over all individual parts
            for part=1:length(simRes(r).t)

                idx = contains(val, simRes(r).t{part}');

                % save location
                loc{end+1,1} = simRes(r).loc(part,:);

                % save corresponding time
                t{end+1,1} = simRes(r).t{part}(idx);

                % save state vector (if given)
                if ~isempty(simRes(r).x)
                    x{end+1,1} = simRes(r).x{part}(idx,:);
                end
                % save output vector (if given)
                if ~isempty(simRes(r).y)
                    y{end+1,1} = simRes(r).y{part}(idx,:);
                end
                % save algebraic vector (if given)
                if ~isempty(simRes(r).a)
                    a{end+1,1} = simRes(r).a{part}(idx,:);
                end
            end

            % append resulting partial simRes object (only if non-empty)
            
            res = [res; simResult(x,t,loc, y, a)];
        end

        
    otherwise
        throw(CORAerror('CORA:wrongValue','second',"'location', 'time'"));

end

% ------------------------------ END OF CODE ------------------------------
