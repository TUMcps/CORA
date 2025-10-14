function res = find(traj,prop,val)
% find - get traj object that satisfy a certain condition
%
% Syntax:
%    res = find(traj,prop,val)
%
% Inputs:
%    traj - trajectory object
%    prop - property for condition ('location', 'time')
%    val - value for property
%
% Outputs:
%    res - all parts of traj trajectories that satisfy the condition
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet/find

% Authors:       Mark Wetzlinger, Tobias Ladner, Laura Luetzow
% Written:       26-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input arguments (only first two)
inputArgsCheck({{traj,'att','trajectory'};
                {prop,'str',{'location', 'time'}}});

% check which property desired
switch prop

    case 'location'
        % get all traj trajectories with a given location number

        % quick check: all trajectories from continuous-time system?
        if isscalar(val) && val == 0 && ...
                all(arrayfun(@(x) all(x.loc==val),traj,'UniformOutput',true))
            res = traj; return
        end

        res = [];
        % loop over all trajectories
        for r=1:length(traj)

            % init properties of trajectory object
            u = []; x = []; y = []; a = []; t = []; loc = [];
            [~,~,~,idzEnd] = extractHits(traj(r),val);
            [~,~,~,idzStart] = extractHits(traj(r),[],val);
            idzStart = idzStart + 1;

            if all(traj(r).loc(:,1) == val)
                % condition is satisfied at the beginning
                idzStart = [1 idzStart];
            end
            if all(traj(r).loc(:,end) == val)
                % condition is satisfied at the end
                idzEnd = [idzEnd traj(r).n_k];
            end
            if isempty(idzStart) && isempty(idzEnd)
                return
            end

            % loop over all individual parts
            for i = 1:length(idzStart)
                iStart = idzStart(i);
                iEnd = idzEnd(i);

                % extract locations
                loc = traj(r).loc(:,iStart:iEnd);

                % extract inputs
                if ~isempty(traj(r).u)
                    u = traj(r).u(:,iStart:iEnd,:);
                end

                % save states (if given)
                if ~isempty(traj(r).x) && size(traj(r).x,2) > 1
                    x = traj(r).x(:,iStart:iEnd,:);
                end
                % save output vector (if given)
                if ~isempty(traj(r).y)
                    y = traj(r).y(:,iStart:iEnd,:);
                end

                % save corresponding time
                if ~isempty(traj(r).t)
                    t = traj(r).t(:,iStart:iEnd);
                end

                % save algebraic vector (if given)
                if ~isempty(traj(r).a)
                    a = traj(r).a(:,iStart:iEnd,:);
                end
                res = [res; trajectory(u,x,y,t,[],loc,a)];
            end
        end

    case 'time'

        % read val
        if ~isa(val,'interval')
            val = interval(val-1e-10,val+1e-10);
        end

        res(length(traj),1) = trajectory();

        % loop over all trajectories
        for r=1:length(traj)

            % init properties of trajectory object
            u = []; x = []; y = []; a = []; t = []; loc = [];

            % loop over all individual parts
            idx = contains(val, traj(r).t);

            % save location
            if ~isempty(traj(r).loc)
                loc = traj(r).loc(:,idx);
            end

            % save corresponding input
            if ~isempty(traj(r).u)
                u = traj(r).u(:,idx,:);
            end
            % save state vector (if given)
            if ~isempty(traj(r).x) && size(traj(r).x, 2) > 1
                x = traj(r).x(:,idx,:);
            end
            % save output vector (if given)
            if ~isempty(traj(r).y)
                y = traj(r).y(:,idx,:);
            end

            % save corresponding time
            if ~isempty(traj(r).t)
                t = traj(r).t(:,idx,:);
            end
            % save algebraic vector (if given)
            if ~isempty(traj(r).a)
                a = traj(r).a(:,idx,:);
            end

            % append resulting partial traj object (only if non-empty)

            res(r) = trajectory(u,x,y,t,[],loc,a);
        end

    otherwise
        throw(CORAerror('CORA:wrongValue','second',"'location', 'time'"));

end

% ------------------------------ END OF CODE ------------------------------
