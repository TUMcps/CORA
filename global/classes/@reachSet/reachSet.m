classdef reachSet
% reachSet - class that stores reachable sets
%
% Syntax:  
%    obj = reachSet(timePoint)
%    obj = reachSet(timePoint,parent)
%    obj = reachSet(timePoint,parent,loc)
%    obj = reachSet(timePoint,timeInt)
%    obj = reachSet(timePoint,timeInt,parent);
%    obj = reachSet(timePoint,timeInt,parent,loc);
%
% Inputs:
%    timePoint - struct with fields .set and .time storing the time point
%                reachable set
%    timeInt - struct with fields .set, .time, and .algebraic (nonlinDAsys)
%              storing the time interval reachable set
%    parent - index of the parent reachable set
%    loc - index of the location (hybrid systems)
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reach

% Author:       Niklas Kochdumper
% Written:      29-May-2020             
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    
    timePoint = [];                     % time point reachable set
    timeInterval = [];                  % time interval reachable set
    parent (1,1) {mustBeInteger} = 0;   % index of parent reachable set
    loc (:,1) {mustBeInteger} = 0;      % index of the location (hybrid)
end
    
methods
    
    % class constructor
    function obj = reachSet(varargin)
        
        % parse input arguments
        if nargin == 1
            obj.timePoint = varargin{1};
        elseif nargin == 2
            obj.timePoint = varargin{1};
            if ~isstruct(varargin{2})
                obj.parent = varargin{2};
            else
                obj.timeInterval = varargin{2};
            end
        elseif nargin == 3
            obj.timePoint = varargin{1};
            if ~isstruct(varargin{2})
                obj.parent = varargin{2};
                obj.loc = varargin{3};
            else
                obj.timeInterval = varargin{2};
                obj.parent = varargin{3};
            end
        elseif nargin == 4
            obj.timePoint = varargin{1};
            obj.timeInterval = varargin{2};
            obj.parent = varargin{3};
            obj.loc = varargin{4};
        else
           error('Wrong number of input arguments for class "reachSet"!'); 
        end  
        
        % check correctness of inputs
        if ~isempty(obj.timeInterval)
           temp = obj.timeInterval;
           if ~isstruct(temp) || ~isfield(temp,'set') || ...
              ~isfield(temp,'time') || ~iscell(temp.set) || ...
              ~iscell(temp.time) || any(size(temp.set) ~= size(temp.time))
               error('Wrong format for input arguments');
           end
        end
        
        if ~isempty(obj.timePoint)
           temp = obj.timePoint;
           if ~isstruct(temp) || ~isfield(temp,'set') || ...
              ~isfield(temp,'time') || ~iscell(temp.set) || ...
              ~iscell(temp.time) || any(size(temp.set) ~= size(temp.time))
               error('Wrong format for input arguments');
           end
        end
    end   
    
    % assign array elements
    function obj = subsasgn(obj, S, value)
        % call built-in function
        obj = builtin('subsasgn', obj, S, value);
    end
    
    % get array entries
    function res = subsref(obj, S)
        % call built-in function
        res = builtin('subsref', obj, S);
    end
end
end

%------------- END OF CODE --------------