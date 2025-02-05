classdef taylorLinSys < handle
% taylorLinSys - helper class for storing auxiliary values used in the
%    reachability analysis of linear continuous-time systems
% note: do not use outside of built-in reachability algorithms
%
% Syntax:
%    obj = taylorLinSys(A)
%
% Inputs:
%    A - state matrix
%
% Outputs:
%    obj - generated taylorLinSys object
%
% Example:
%    A = [-2 0; 1 -3];
%    sys = taylorLinSys(A)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys

% Authors:       Mark Wetzlinger
% Written:       03-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public) %(Access = private)
    A;                  % system matrix (square)
    A_abs;              % absolute value of system matrix
    Ainv;               % inverse of system matrix
    timeStep;           % array of time steps that have been used

    % fields depending on timeStep
    eAdt;               % propagation matrix e^(A*dt)
    E;                  % remainder of exponential matrix
    F;                  % correction matrix for the state
    G;                  % correction matrix for the input
    dtoverfac;          % factor dt^i/i!

    % fields depending on truncationOrder (directly as index usable)
    % note: zeroth element is always identity matrix, so index 1 stores i=1
    Apower;             % cell-array of A^i: A, A^2, A^3, ...
    Apower_abs;         % cell-array of |A|^i: |A|, |A|^2, |A|^3, ...
    Apos;               % cell-array of indices of positive values of A^i
    Aneg;               % cell-array of indices of negative values of A^i

    % fields depending on timeStep (used as first index) and
    % truncationOrder (used as second index) -- this helps in keeping
    % rounding errors small, as powers of A converge to Inf or 0 and dt^i/i!
    % converges to 0
    Apower_dt_fact;     % cell-array of (A*dt)^i/i!: A*dt, (A*dt)^2/2, (A*dt)^3/3!, ...
    Apower_abs_dt_fact; % cell-array of (|A|*dt)^i/i!: |A|*dt, (|A|*dt)^2/2, (|A|*dt)^3/3!, ...
end

methods
    % constructor
    function obj = taylorLinSys(A)
        % brief check
        if ~isnumeric(A) || numel(size(A)) > 2 || size(A,1) ~= size(A,2) ...
                || ~all(isfinite(A),'all') || ~all(isreal(A),'all')
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'A must be numeric, finite, 2D, square, and real-valued'));
        end

        obj.A = A;
        obj.A_abs = abs(A);
        % don't compute inverse because A might be singular... (only if
        % it is explicitly requested)

        % init first entry (often outside of loops)
        obj.Apower{1} = A;
        obj.Apower_abs{1} = abs(A);
        obj.Apos{1} = A;
        obj.Apos{1}(A < 0) = 0;
        obj.Aneg{1} = A;
        obj.Aneg{1}(A > 0) = 0;
        % init type
        obj.Apower_dt_fact = {};
        obj.Apower_abs_dt_fact = {};
    end

    function val = computeField(obj,name,options)
        % options with fields 'timeStep', 'ithpower'
        switch name
            case 'Ainv'
                val = compute_Ainv(obj);
            case 'eAdt'
                val = compute_eAdt(obj,options.timeStep);
            case 'Apower'
                val = compute_Apower(obj,options.ithpower);
            case 'Apower_abs'
                val = compute_Apower_abs(obj,options.ithpower);
            case 'Apower_dt_fact'
                val = compute_Apower_dt_fact(obj,options.timeStep,options.ithpower);
            case 'Apower_abs_dt_fact'
                val = compute_Apower_abs_dt_fact(obj,options.timeStep,options.ithpower);
            case 'Apos'
                % note: positive indices not affected by * dt/factorial
                val = compute_Apos(obj,options.ithpower);
            case 'Aneg'
                % note: negative indices not affected by * dt/factorial
                val = compute_Aneg(obj,options.ithpower);
            case 'dtoverfac'
                val = compute_dtoverfac(obj,options.timeStep,options.ithpower);
            otherwise
                throw(CORAerror('CORA:wrongValue','second',...
                        ['has to be "Ainv", "eAdt", "Apower", "Apower_abs", ' ...
                        '"Apower_dt_fact", "Apower_abs_dt_fact", "Apos", '...
                        '"Aneg", or "dtoverfac"']));
        end
    end

    function Apower_i = compute_Apower(obj,ithpower)
        % check if already there
        if length(obj.Apower) >= ithpower
            Apower_i = readField(obj,'Apower',ithpower);
            return
        end

        % compute and save new value
        Apower_mm = readField(obj,'Apower',ithpower-1);
        if isempty(Apower_mm)
            % recursive call (idx = 1 always computed, see constructor)
            Apower_mm = compute_Apower(obj,ithpower-1);
        end
        Apower_i = Apower_mm * obj.A;
        obj.Apower{ithpower} = Apower_i;
    end
    function Apower_abs_i = compute_Apower_abs(obj,ithpower)
        % check if already there
        if length(obj.Apower_abs) >= ithpower
            Apower_abs_i = readField(obj,'Apower_abs',ithpower);
            return
        end

        % compute and save new value
        Apower_abs_mm = readField(obj,'Apower_abs',ithpower-1);
        if isempty(Apower_abs_mm)
            % recursive call (idx = 1 always computed, see constructor)
            Apower_abs_mm = compute_Apower_abs(obj,ithpower-1);
        end
        Apower_abs_i = Apower_abs_mm * obj.A_abs;
        obj.Apower_abs{ithpower} = Apower_abs_i;
    end
    function Apower_dt_fact_val = compute_Apower_dt_fact(obj,timeStep,ithpower)
        % check if time step has been used before
        timeStepIdx = getIndexForTimeStep(obj,timeStep);
        if timeStepIdx == -1
            makeNewTimeStep(obj,timeStep);
            timeStepIdx = length(obj.timeStep);
            Apower_dt_fact_cell = [];
        else
            Apower_dt_fact_cell = obj.Apower_dt_fact{timeStepIdx};
        end
        % check if already there
        if ~isempty(Apower_dt_fact_cell) && numel(Apower_dt_fact_cell) >= ithpower
            Apower_dt_fact_val = readField(obj,'Apower_dt_fact',timeStepIdx,ithpower);
            return
        end

        % compute and save new power
        if isempty(Apower_dt_fact_cell) && ithpower == 1
            % initialize cell array for given time step
            Apower_dt_fact_val = obj.A*timeStep;
            obj.Apower_dt_fact{timeStepIdx}{ithpower} = Apower_dt_fact_val;
        elseif numel(Apower_dt_fact_cell) < ithpower-1
            % recursive call until ithpower-1 is computed
            Apower_dt_fact_val = compute_Apower_dt_fact(obj,timeStep,ithpower-1);
            Apower_dt_fact_val = Apower_dt_fact_val * obj.A * timeStep / ithpower;
            obj.Apower_dt_fact{timeStepIdx}{ithpower} = Apower_dt_fact_val;
        else
            % last one at ithpower-1 is computed
            Apower_dt_fact_val = Apower_dt_fact_cell{ithpower-1} * obj.A * timeStep / ithpower;
            obj.Apower_dt_fact{timeStepIdx}{ithpower} = Apower_dt_fact_val;
        end
    end
    function Apower_abs_dt_fact_val = compute_Apower_abs_dt_fact(obj,timeStep,ithpower)
        % check if time step has been used before
        timeStepIdx = getIndexForTimeStep(obj,timeStep);
        if timeStepIdx == -1
            makeNewTimeStep(obj,timeStep);
            timeStepIdx = length(obj.timeStep);
            Apower_abs_dt_fact_cell = [];
        else
            Apower_abs_dt_fact_cell = obj.Apower_abs_dt_fact{timeStepIdx};
        end
        % check if already there
        if ~isempty(Apower_abs_dt_fact_cell) && numel(Apower_abs_dt_fact_cell) >= ithpower
            Apower_abs_dt_fact_val = readField(obj,'Apower_abs_dt_fact',timeStepIdx,ithpower);
            return
        end

        % compute and save new power
        if isempty(Apower_abs_dt_fact_cell) && ithpower == 1
            % initialize cell array for given time step
            Apower_abs_dt_fact_val = obj.A_abs*timeStep;
            obj.Apower_abs_dt_fact{timeStepIdx}{ithpower} = Apower_abs_dt_fact_val;
        elseif numel(Apower_abs_dt_fact_cell) < ithpower-1
            % recursive call until ithpower-1 is computed
            Apower_abs_dt_fact_val = compute_Apower_abs_dt_fact(obj,timeStep,ithpower-1);
            Apower_abs_dt_fact_val = Apower_abs_dt_fact_val * obj.A_abs * timeStep / ithpower;
            obj.Apower_abs_dt_fact{timeStepIdx}{ithpower} = Apower_abs_dt_fact_val;
        else
            % last one at ithpower-1 is computed
            Apower_abs_dt_fact_val = Apower_abs_dt_fact_cell{ithpower-1} * obj.A_abs * timeStep / ithpower;
            obj.Apower_abs_dt_fact{timeStepIdx}{ithpower} = Apower_abs_dt_fact_val;
        end
    end
    function Apos = compute_Apos(obj,ithpower)
        % check if already there
        if length(obj.Apos) >= ithpower
            Apos = readField(obj,'Apos',ithpower);
            return
        end

        % compute/read out Apower_eta
        Apower_i = compute_Apower(obj,ithpower);
        % get positive indices
        Apos = Apower_i;
        Apos(Apower_i < 0) = 0;
        % save
        obj.Apos{ithpower} = Apos;
    end
    function Aneg = compute_Aneg(obj,ithpower)
        % check if already there
        if length(obj.Aneg) >= ithpower
            Aneg = readField(obj,'Aneg',ithpower);
            return
        end

        % compute/read out Apower_eta
        Apower_i = compute_Apower(obj,ithpower);
        % get negative indices
        Aneg = Apower_i;
        Aneg(Apower_i > 0) = 0;
        % save
        obj.Aneg{ithpower} = Aneg;
    end
    function dtoverfac = compute_dtoverfac(obj,timeStep,ithpower)
        % check if given time step has been used before
        idx = getIndexForTimeStep(obj,timeStep);
        if idx == -1
            makeNewTimeStep(obj,timeStep);
            idx = length(obj.timeStep);
            dtoverfac_mm = [];
        else
            dtoverfac_mm = obj.dtoverfac{idx};
        end
        
        % check if already there
        if ~isempty(dtoverfac_mm) && length(dtoverfac_mm) >= ithpower
            dtoverfac = obj.dtoverfac{idx}(ithpower);
            return
        end

        % compute and save new power
        if isempty(dtoverfac_mm) && ithpower == 1
            dtoverfac = timeStep;
            obj.dtoverfac{idx}(ithpower) = dtoverfac;
        elseif length(dtoverfac_mm) < ithpower-1
            % recursive call (idx = 1 always computed, see constructor)
            dtoverfac = compute_dtoverfac(obj,timeStep,ithpower-1);
            dtoverfac = dtoverfac * timeStep / ithpower;
            obj.dtoverfac{idx}(ithpower) = dtoverfac;
        else
            dtoverfac = dtoverfac_mm(end) * timeStep / ithpower;
            obj.dtoverfac{idx}(ithpower) = dtoverfac;
        end
    end
    function eAdt = compute_eAdt(obj,timeStep)
        % compute eAdt
        eAdt = readFieldForTimeStep(obj,'eAdt',timeStep);
        if isempty(eAdt)
            eAdt = expm(obj.A*timeStep);
            insertFieldTimeStep(obj,'eAdt',eAdt,timeStep);
        end
    end
    function val = compute_Ainv(obj)
        % compute Ainv
        if ~isempty(obj.Ainv)
            val = obj.Ainv;
            return;
        end

        if rank(obj.A) < length(obj.A)
            val = [];
        else
            val = inv(obj.A);
        end
        obj.Ainv = val;
    end

    % read out (only if already there)
    function val = readFieldForTimeStep(obj,field,timeStep)
        idx = getIndexForTimeStep(obj,timeStep);
        val = readField(obj,field,idx);
    end
    function val = readField(obj,field,idx1,idx2)
        % idx1 relates to timeStep, idx2 to ithpower
        val = [];
        if idx1 ~= -1
            try
                switch field
                    case 'E'
                        val = obj.E{idx1};
                    case 'F'
                        val = obj.F{idx1};
                    case 'G'
                        val = obj.G{idx1};
                    case 'eAdt'
                        val = obj.eAdt{idx1};
                        % A
                    case 'Apower'
                        val = obj.Apower{idx1};
                    case 'Apower_abs'
                        val = obj.Apower_abs{idx1};
                    case 'Apower_dt_fact'
                        val = obj.Apower_dt_fact{idx1}{idx2};
                    case 'Apower_abs_dt_fact'
                        val = obj.Apower_abs_dt_fact{idx1}{idx2};
                    case 'Apos'
                        val = obj.Apos{idx1};
                    case 'Aneg'
                        val = obj.Aneg{idx1};
                    case 'dtoverfac'
                        val = obj.dtoverfac{idx1};
                    otherwise % throw error
                        throw(CORAerror('CORA:wrongValue','third',...
                            ['has to be "E", "F", "G", "Apower", "Apower_abs", ' ...
                            '"Apos", "Aneg", or "dtoverfac"']));
                end
            catch
                % likely index out of range
                val = [];
            end
        end
    end

    % important helper indexing function
    function idx = getIndexForTimeStep(obj,timeStep)
        % check if a given time step size has already been used, return the
        % corresponding index in the list; (-1 if not there)
        idx = -1;
        idxLogical = withinTol(obj.timeStep,timeStep,1e-10);
        if any(idxLogical)
            idx = find(idxLogical);
        end
    end
    function makeNewTimeStep(obj,timeStep)
        % append time step to list of time steps
        obj.timeStep(end+1,1) = timeStep;
        % increase length of all fields
        obj.eAdt{end+1,1} = [];
        obj.Apower_dt_fact{end+1,1} = [];
        obj.Apower_abs_dt_fact{end+1,1} = [];
        obj.E{end+1,1} = [];
        obj.F{end+1,1} = [];
        obj.G{end+1,1} = [];
        obj.dtoverfac{end+1,1} = [];
    end

    function insertFieldTimeStep(obj,field,val,timeStep)
        idx = getIndexForTimeStep(obj,timeStep);
        if idx == -1
            makeNewTimeStep(obj,timeStep);
            idx = length(obj.timeStep);
        end
        insertField(obj,field,val,idx);
    end
    function insertField(obj,field,val,idx)
        % check if already there?
        switch field
            case 'E'
               obj.E{idx} = val;
            case 'F'
                obj.F{idx} = val;
            case 'G'
                obj.G{idx} = val;
            case 'eAdt'
                obj.eAdt{idx} = val;
            otherwise
                throw(CORAerror('CORA:wrongValue','third',...
                        'has to be "E", "F", "G", or "eAdt"'));
        end
    end

end

end

% ------------------------------ END OF CODE ------------------------------
