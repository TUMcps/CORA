classdef simResult
% simResult - class that stores simulation results
%
% Syntax:
%    simRes = simResult(x,t)
%    simRes = simResult(x,t,loc)
%    simRes = simResult(x,t,{},y)
%    simRes = simResult(x,t,{},y,a)
%
% Inputs:
%    x - cell-array storing the states of the simulated trajectories, where
%        each trajectory is a matrix of dimension [N,n]
%    t - cell-array storing the time points for the simulated trajectories
%        as vectors of dimension [N,1]
%    loc - double-array storing the locations for the simulated trajectories
%    y - cell-array storing the output of the simulated trajectories, where
%        each trajectory is a matrix of dimension [N,o]
%    a - cell-array storing the algebraic variables of the simulated
%        trajectories (only for nonlinDASys), where each trajectory is a
%        matrix of dimension [N,p]
%    (where N ... number of simulated trajectories, n ... state dimension,
%     o ... output dimension, p ... number of algebraic variables)
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reachSet, simulateRandom

% Authors:       Niklas Kochdumper
% Written:       29-May-2020
% Last update:   16-November-2021 (MW, add property .y)
%                02-June-2022 (MW, add property .a)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    
    x = {};             % states of the simulated trajectories
    y = {};             % outputs of the simulated trajectories
    a = {};             % algebraic variables of the sim. trajectories
    t = {};             % time of the simulated trajectories
    loc = 0;            % index of the locations (0 for contDynamics)
    
end


methods
    
    % class constructor
    function simRes = simResult(varargin)
        
        % parse input arguments
        if nargin == 0
            % empty object
            return
        elseif nargin == 1 && isa(varargin{1},'simResult')
            % copy constructor
            simRes = varargin{1};
        elseif nargin == 1
            throw(CORAerror('CORA:notEnoughInputArgs',2));
        elseif nargin == 2
            simRes.x = varargin{1};
            simRes.t = varargin{2};
        elseif nargin == 3
            simRes.x = varargin{1};
            simRes.t = varargin{2};
            simRes.loc = varargin{3};
        elseif nargin == 4
            simRes.x = varargin{1};
            simRes.t = varargin{2};
            simRes.loc = varargin{3};
            simRes.y = varargin{4};
        elseif nargin == 5
            simRes.x = varargin{1};
            simRes.t = varargin{2};
            simRes.loc = varargin{3};
            simRes.y = varargin{4};
            simRes.a = varargin{5};
        else
        	throw(CORAerror('CORA:tooManyInputArgs',5));
        end


        % check correct format of simulated trajectories
        if CHECKS_ENABLED
            % - time has to be non-empty, the entries have to be nonnegative
            %   and monotonically increasing
            if isempty(simRes.t) || ~iscell(simRes.t) || ...
                    any(cellfun(@(x) size(x,2),simRes.t,'UniformOutput',true) ~= 1) || ...
                    any(cellfun(@(x) any(x < 0),simRes.t,'UniformOutput',true)) || ...
                    any(cellfun(@(x) any(diff(x < 0)),simRes.t,'UniformOutput',true))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Property .t has the wrong format.'));
            end
            tVeclength = cellfun(@(x) length(x),simRes.t,'UniformOutput',true);
            % check which trajectories are given
            props = {'x','y','a'};
            isTraj = [~isempty(simRes.x), ~isempty(simRes.y), ~isempty(simRes.a)];
    
            % all trajectories have to have same number of runs
            if nnz(isTraj) > 1
                lengths = [length(simRes.x) length(simRes.y) length(simRes.a)];
                if any(diff(lengths(isTraj)) ~= 0)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'All trajectories have to have the same number of runs.'));
                end
            end
    
            % number of steps of each run has to equal length of time vector
            for i=1:length(props)
                proplength = cellfun(@(x) size(x,1),simRes.(props{i}),'UniformOutput',true);
                if ~isempty(proplength) && any(proplength ~= tVeclength)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'The length of each run has to match the length of the time vector.'));
                end
            end
    
            % number of columns has to be the same for each group of trajectories
            for i=1:length(props)
                if isTraj(i) && length(simRes.(props{i})) > 1
                    cols = cellfun(@(x) size(x,2),simRes.(props{i}),'UniformOutput',true);
                    if any(diff(cols) > 0)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'The number of columns in .x/.y,/.a has to be the same in each run'));
                    end
                end
            end
        end

    end
end

end

% ------------------------------ END OF CODE ------------------------------
