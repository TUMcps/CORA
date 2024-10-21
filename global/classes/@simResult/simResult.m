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
%                29-June-2024 (TL, allow empty .t if .x is empty)
% Last revision: 15-October-2024 (MW, restructure)

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

        % 0. avoid empty instantiation
        if nargin == 0
            % empty object
            return
        end
        assertNarginConstructor(1:5,nargin);

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'simResult')
            % copy constructor
            simRes = varargin{1};
        end
        
        % 2. parse input arguments: varargin -> vars
        [x,t,loc,y,a] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(x,t,loc,y,a,nargin);

        % 4. assign properties
        simRes.x = x;
        simRes.t = t;
        simRes.loc = loc;
        simRes.y = y;
        simRes.a = a;

    end
end

end


% Auxiliary functions -----------------------------------------------------

function [x,t,loc,y,a] = aux_parseInputArgs(varargin)

    x = {}; t = {}; loc = 0; y = {}; a = {};
    [x,t,loc,y,a] = setDefaultValues({x,t,loc,y,a},varargin);

end

function aux_checkInputArgs(x,t,loc,y,a,n_in)

% check correct format of simulated trajectories
if CHECKS_ENABLED && n_in > 0
    % - time has to be non-empty, the entries have to be nonnegative
    %   and monotonically increasing
    if (isempty(t) && ~isempty(x)) || ~iscell(t) || ...
            any(cellfun(@(x_) size(x_,2),t,'UniformOutput',true) ~= 1) || ...
            any(cellfun(@(x_) any(x_ < 0),t,'UniformOutput',true)) || ...
            any(cellfun(@(x_) any(diff(x_ < 0)),t,'UniformOutput',true))
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Property .t has the wrong format.'));
    end
    tVeclength = cellfun(@(x_) length(x_),t,'UniformOutput',true);
    % check which trajectories are given
    props = {x,y,a};
    isTraj = [~isempty(x), ~isempty(y), ~isempty(a)];

    % all trajectories have to have same number of runs
    if nnz(isTraj) > 1
        lengths = [length(x) length(y) length(a)];
        if any(diff(lengths(isTraj)) ~= 0)
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'All trajectories have to have the same number of runs.'));
        end
    end

    % number of steps of each run has to equal length of time vector
    for i=1:length(props)
        proplength = cellfun(@(x_) size(x_,1),props{i},'UniformOutput',true);
        if ~isempty(proplength) && any(proplength ~= tVeclength)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The length of each run has to match the length of the time vector.'));
        end
    end

    % number of columns has to be the same for each group of trajectories
    for i=1:length(props)
        if isTraj(i) && length(props{i}) > 1
            cols = cellfun(@(x_) size(x_,2),props{i},'UniformOutput',true);
            if any(diff(cols) > 0)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'The number of columns in .x/.y,/.a has to be the same in each run'));
            end
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
