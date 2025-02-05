classdef reachSet
% reachSet - class that stores reachable sets
%
% Syntax:
%    R = reachSet(timePoint)
%    R = reachSet(timePoint,parent)
%    R = reachSet(timePoint,parent,loc)
%    R = reachSet(timePoint,timeInt)
%    R = reachSet(timePoint,timeInt,parent)
%    R = reachSet(timePoint,timeInt,parent,loc)
%
% Inputs:
%    timePoint - struct with fields .set, .time, and .error (linearSys if
%                adaptive algorithm called) storing the time point
%                reachable or output set
%    timeInt - struct with fields .set, .time, .error (linearSys if
%              adaptive algorithm called), and .algebraic (nonlinDASys)
%              storing the time interval reachable or output set
%    parent - index of the parent reachable set
%    loc - index of the location (hybrid systems)
%
% Outputs:
%    R - generated reachSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reach

% Authors:       Niklas Kochdumper
% Written:       29-May-2020             
% Last update:   ---
% Last revision: 15-October-2024 (MW, restructure)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    
    % time-point reachable set
    timePoint = [];

    % time-interval reachable set
    timeInterval = [];

    % index of parent reachable set
    parent (1,1) {mustBeInteger,mustBeNonnegative} = 0;

    % index of the location (hybrid)
    loc (:,1) {mustBeInteger,mustBeNonnegative} = 0;
end
    
methods
    
    % class constructor
    function R = reachSet(varargin)
        
        % 1. empty instantiation
        if nargin == 0
            return
        end
        assertNarginConstructor(1:4,nargin);

        % 2. copy constructor
        if nargin == 1
            if isa(varargin{1},'reachSet')
                R = varargin{1};
                return
            end
        end
        
        % 2. parse input arguments: varargin -> vars
        [timePoint,timeInterval,parent,loc] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(timePoint,timeInterval,parent,loc,nargin);

        % 4. assign properties
        R.timePoint = timePoint;
        R.timeInterval = timeInterval;
        R.parent = parent;
        R.loc = loc;

    end

    % get initial set
    function R0 = R0(obj)
        if isempty(obj(1).timePoint)
            n = 3; % ensures successful plotting
            R0 = initialSet(emptySet(n));
        else
            R0 = initialSet(obj(1).timePoint(1).set{1});
        end
    end

    % for functionSignatures.json
    han = plot(obj, dims, varargin)
    han = plotOverTime(obj, dim, varargin)
end

methods (Static = true)
    % instantiate reachSet object from structs
    R = initReachSet(timePoint,timeInt)
end

end


% Auxiliary functions -----------------------------------------------------

function [timePoint,timeInterval,parent,loc] = aux_parseInputArgs(varargin)

    % default values
    timePoint = []; timeInterval = [];

    % timePoint given?
    if numel(varargin) >= 1
        timePoint = varargin{1};
    end

    % timeInterval given?
    if numel(varargin) >= 2 && ... % check condition for timeInterval
        (isstruct(varargin{2}) || (isnumeric(varargin{2}) && isempty(varargin{2})))
        % read timeInterval and remove from varargin
        timeInterval = varargin{2};
        varargin = varargin([1,3:end]);
    end

    % read parent and loc
    [parent,loc] = setDefaultValues({0,0},varargin(2:end));

end

function aux_checkInputArgs(timePoint,timeInterval,parent,loc,n_in)

% check correct format of reachable sets
if CHECKS_ENABLED && n_in > 0

    % same checks for both sets
    sets = {timeInterval,timePoint};

    for i=1:numel(sets)
        if isempty(sets{i})
            continue
        end

        Rset = sets{i};
        if ~isstruct(Rset) ... % has to be a struct
                || ~isfield(Rset,'set') || ... % has to have .set
                ~isfield(Rset,'time') || ... % has to have .time
                (~isempty(Rset.set) && ~iscell(Rset.set)) || ... % .set have to be cells
                (~isempty(Rset.time) && ~iscell(Rset.time)) || ... % .time have to be cells
                any(size(Rset.set) ~= size(Rset.time)) % .set and .time have to be of equal length
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Fields .set and/or .time have the wrong format.'));
        end
        % optional fields: .error, .algebraic
        if isfield(Rset,'error') && ...
                ( ~isnumeric(Rset.error) || ... % has to be numeric
                  any(size(Rset.error) ~= size(Rset.time)) ) % correct length
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Field .error have the wrong format.'));
        elseif isfield(Rset,'algebraic') && ...
                ( ~iscell(Rset.algebraic) || ... % has to be cell-array
                  any(size(Rset.algebraic) ~= size(Rset.time)) ) % correct length
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Field .algebraic have the wrong format.'));
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
