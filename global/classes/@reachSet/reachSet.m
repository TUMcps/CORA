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
% Last revision: ---

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
        
        % parse input arguments
        if nargin == 0
            return
            % empty object
        elseif nargin == 1
            if isa(varargin{1},'reachSet')
                % copy constructor
                R = varargin{1};
            else
                R.timePoint = varargin{1};
            end
        elseif nargin == 2
            R.timePoint = varargin{1};
            if isstruct(varargin{2}) || ...
                    (isnumeric(varargin{2}) && isempty(varargin{2}))
                R.timeInterval = varargin{2};
            else
                R.parent = varargin{2};
            end
        elseif nargin == 3
            R.timePoint = varargin{1};
            if isstruct(varargin{2}) || ...
                    (isnumeric(varargin{2}) && isempty(varargin{2}))
                R.timeInterval = varargin{2};
                R.parent = varargin{3};
            else
                R.parent = varargin{2};
                R.loc = varargin{3};
            end
        elseif nargin == 4
            R.timePoint = varargin{1};
            R.timeInterval = varargin{2};
            R.parent = varargin{3};
            R.loc = varargin{4};
        else
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end
        
        % check correct format of reachable sets
        if CHECKS_ENABLED
            sets = {'timeInterval','timePoint'};
            for i=1:length(sets)
                if ~isempty(R.(sets{i}))
                    temp = R.(sets{i});
                    if ~isstruct(temp) ... % has to be a struct
                            || ~isfield(temp,'set') || ... % has to have .set
                            ~isfield(temp,'time') || ... % has to have .time
                            (~isempty(temp.set) && ~iscell(temp.set)) || ... % .set have to be cells
                            (~isempty(temp.time) && ~iscell(temp.time)) || ... % .time have to be cells
                            any(size(temp.set) ~= size(temp.time)) % .set and .time have to be of equal length
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Fields .set and/or .time have the wrong format.'));
                    end
                    % optional fields: .error, .algebraic
                    if isfield(temp,'error') && ...
                            ( ~isnumeric(temp.error) || ... % has to be numeric
                              any(size(temp.error) ~= size(temp.time)) ) % correct length
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Field .error have the wrong format.'));
                    elseif isfield(temp,'algebraic') && ...
                            ( ~iscell(temp.algebraic) || ... % has to be cell-array
                              any(size(temp.algebraic) ~= size(temp.time)) ) % correct length
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Field .algebraic have the wrong format.'));
                    end
                end
            end
        end
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

% ------------------------------ END OF CODE ------------------------------
