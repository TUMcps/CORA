classdef specification
% specfication - class that stores spefications
%
% Syntax:  
%    obj = specification(set)
%    obj = specification(list)
%    obj = specification(set,type)
%    obj = specification(list,type)
%    obj = specification(func,'custom')
%
% Inputs:
%    set - contSet object that defines the specification
%    value - string that defines the type of spefication. Supported types
%            are 'unsafeSet' (default), 'safeSet', 'custom', or 'invariant' 
%    list - cell-array storing with contSet objects for multiple parallel
%           specifications
%    func - function handle to a user provided specification check function
%
% Outputs:
%    obj - generated object
%
% Example:
%    h = halfspace([1,2],0);
%    spec = specification(h,'unsafeSet');
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
    
    % contSet object that corresponds to the specification
    set = [];
    
    % type of specification
    type (1,:) char {mustBeMember(type, ...
         {'unsafeSet','safeSet','invariant','custom'})} = 'unsafeSet';
end
    
methods
    
    % class constructor
    function obj = specification(varargin)
        
        % parse input arguments
        if nargin == 1
            if iscell(varargin{1})
                obj = repelem(obj,length(varargin{1}),1);
                for i = 1:length(varargin{1})
                   obj(i,1).set = varargin{1}{i};
                end
            else
                obj.set = varargin{1};
            end
        elseif nargin == 2
            if iscell(varargin{1})
                obj = repelem(obj,length(varargin{1}),1);
                for i = 1:length(varargin{1})
                   obj(i,1).set = varargin{1}{i};
                   obj(i,1).type = varargin{2};
                end
            else
                obj.set = varargin{1};
                obj.type = varargin{2};
            end
        else
           error('Wrong number of input arguments for class "reachSet"!'); 
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