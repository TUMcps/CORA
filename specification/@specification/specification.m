classdef specification
% specification - class that stores specifications
%
% Syntax:
%    obj = specification()
%
%    % single set
%    obj = specification(set)
%    obj = specification(set,location)
%    obj = specification(set,type)
%    obj = specification(set,type,location)
%    obj = specification(set,type,time)
%    obj = specification(set,type,location,time)
%    obj = specification(set,type,time,location)
%
%    % list of sets
%    obj = specification(list)
%    obj = specification(list,location)
%    obj = specification(list,type)
%    obj = specification(list,type,location)
%    obj = specification(list,type,time)
%    obj = specification(list,type,location,time)
%    obj = specification(list,type,time,location)
%
%    % special case: function handle
%    obj = specification(func)
%    obj = specification(func,'custom')
%    obj = specification(func,'custom',location)
%    obj = specification(func,'custom',time)
%    obj = specification(func,'custom',location,time)
%    obj = specification(func,'custom',time,location)
%
%    % special case: stl formula
%    obj = specification(eq)
%    obj = specification(eq,'logic')
%
% Inputs:
%    set - contSet object that defines the specification
%    list - cell-array storing with contSet objects for multiple parallel
%           specifications
%    type - string that defines the type of spefication:
%               - 'unsafeSet' (default)
%               - 'safeSet'
%               - 'invariant' 
%               - 'logic'
%               - 'custom'
%    time - interval defining when the specification is active
%    eq - temporal logic formula (class stl)
%    func - function handle to a user-provided specification check function
%    location - activity of specification in which locations of a HA/pHA
%
% Outputs:
%    obj - generated specification object
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

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       29-May-2020
% Last update:   27-November-2022 (MW, add location property and checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    
    % contSet object that corresponds to the specification
    set = [];
    
    % time interval in which the specification is active
    time interval = []; 
    
    % type of specification
    type (1,:) char {mustBeMember(type, ...
      {'unsafeSet','safeSet','invariant','custom','logic'})} = 'unsafeSet';

    % location where the specification is active (only hybrid systems)
    location = [];
end

methods
    
    % class constructor
    function obj = specification(varargin)
        
        % parse input arguments
        if nargin > 4
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end

        if nargin == 0
            % empty object
            return
        elseif nargin >= 1
            if isa(varargin{1},'specification')
                % copy constructor
                obj = varargin{1}; return
            % first input argument: func, eq, set, list
            elseif isa(varargin{1},'function_handle')
                % syntax: obj = specification(func)
                obj.set = varargin{1};
                obj.type = 'custom';
                isFunHan = true;
                isStl = false;
            elseif isa(varargin{1},'stl')
                % syntax: obj = specification(eq)
                obj.set = varargin{1};
                obj.type = 'logic';
                isFunHan = false;
                isStl = true;
            elseif isa(varargin{1},'contSet')
                % syntax: obj = specification(set)
                obj.set = varargin{1};
                isFunHan = false;
                isStl = false;
            elseif iscell(varargin{1})
                % syntax: obj = specification(list)
                % ensure that list is made of contSet object
                isContSet = all(cellfun(@(x) isa(x,'contSet'),varargin{1},'UniformOutput',true));
                isStl = all(cellfun(@(x) isa(x,'stl'),varargin{1},'UniformOutput',true));
                isFunHan = all(cellfun(@(x) isa(x,'function_handle'),varargin{1},'UniformOutput',true));
                if ~isContSet && ~isStl && ~isFunHan
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        ['Input argument ''list'' has to be a cell-array '...
                        'of contSet/stl/function_handle objects.']));
                end
                % init nx1 specification objects
                obj = repelem(obj,length(varargin{1}),1);
                for i = 1:length(varargin{1})
                    obj(i,1).set = varargin{1}{i};
                    % set correct type
                    if isStl
                        obj(i,1).type = 'logic';
                    elseif isFunHan
                        obj(i,1).type = 'custom';
                    end
                end
            end
        end

        % ...if list was given, we already have an nx1 specification object

        if nargin >= 2
            % second input argument: type, location
            % syntax: obj = specification(set,location)
            %         obj = specification(set,type)
            %         obj = specification(list,location)
            %         obj = specification(list,type)
            %         obj = specification(eq,'logic')
            %         obj = specification(func,'custom')
            [type,loc] = aux_readOutInputArg(varargin{2},2);
            
            if ~isempty(type)             % type is given
                % ensure correct types for func, eq
                if isFunHan && ~strcmp(type,'custom')
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        ['If the specification is defined using a function handle, '...
                        'the property ''type'' must be ''custom''.']));
                elseif isStl && ~strcmp(type,'logic')
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        ['If the specification is defined using an stl formula, '...
                        'the property ''type'' must be ''logic''.']));
                end
                % assign value (checked via property validation)
                for i=1:length(obj)
                    obj(i,1).type = type;
                end

            elseif ~isempty(loc)          % location is given
                % not supported for stl formulae
                if isStl
                    throw(CORAerror('CORA:notSupported',...
                        'Specifications using stl formulae not supported for hybrid systems.'));
                end

                % check that format is correct
                aux_checkLocation(loc);

                % assign value
                for i=1:length(obj)
                    obj(i,1).location = loc;
                end
            end
        end

        if nargin >= 3
            % third input argument: location, time
            % syntax: obj = specification(set,type,location)
            %         obj = specification(set,type,time)
            %         obj = specification(list,type,location)
            %         obj = specification(list,type,time)
            %         obj = specification(func,'custom',location)
            %         obj = specification(func,'custom',time)
            [~,loc,time] = aux_readOutInputArg(varargin{3},3);

            % neither time nor location supported for stl formulae
            if isStl
                throw(CORAerror('CORA:notSupported',...
                    ['Specifications using stl formulae not supported ' ...
                    'for hybrid systems or combined with additional ''time'' input.']));
            end

            if ~isempty(loc)                           % location is given
                % check that format is correct
                aux_checkLocation(loc);

                % assign value
                for i=1:length(obj)
                    obj(i,1).location = loc;
                end
            elseif ~representsa_(time,'emptySet',eps)  % time is given
                % assign value
                for i=1:length(obj)
                    obj(i,1).time = time;
                end
            end
        end

        if nargin == 4
            % fourth input argument: location, time
            % syntax: obj = specification(set,type,location,time)
            %         obj = specification(set,type,time,location)
            %         obj = specification(list,type,location,time)
            %         obj = specification(list,type,time,location)
            %         obj = specification(func,'costum',location,time)
            %         obj = specification(func,'costum',time,location)
            [~,loc,time] = aux_readOutInputArg(varargin{4},4);

            if ~isempty(loc)          % location is given
                % check that format is correct
                aux_checkLocation(loc);

                % assign value
                for i=1:length(obj)
                    obj(i,1).location = loc;
                end
            elseif ~isempty(time)     % time is given
                % assign value
                for i=1:length(obj)
                    obj(i,1).time = time;
                end
            end
        end
    end   

end
end


% Auxiliary functions -----------------------------------------------------

function [type,loc,time] = aux_readOutInputArg(argIn,idx)
% checks whether given input argument is
% - type (has to be char)
% - location (has to be a cell-array)
% - time (has to be an interval object)
% ...otherwise an error is thrown

% init as empty
type = []; loc = []; time = [];

if ischar(argIn)
    type = argIn;
elseif isnumeric(argIn) || iscell(argIn)
    % numeric for HA, cell for pHA
    loc = argIn;
elseif isa(argIn,'interval')
    time = argIn;
else
    if idx == 2
        throw(CORAerror('CORA:wrongInputInConstructor',...
            ['The second input argument has to be either a char-array (type) '...
            'or a double-array/cell-array (location).']));
    elseif idx == 3
        throw(CORAerror('CORA:wrongInputInConstructor',...
            ['The third input argument has to be either a double-array/cell-array (location) '...
            'or an interval object (time).']));
    elseif idx == 4
        throw(CORAerror('CORA:wrongInputInConstructor',...
            ['The fourth input argument has to be either a double-array/cell-array (location) '...
            'or an interval object (time).']));
    end
end

end

function aux_checkLocation(loc)
% checks if the location property has the correct format
% - HA: double array or empty, e.g.,
%        [1,2] = active in locations 1 and 2
% - pHA: cell-array of doubles, e.g.,
%        { [1,2], [1,3] }
%        = active in locations 1 and 2 of subcomponent 1,
%             and in locations 1 and 3 of subcomponent 2

% check whether HA or pHA given
if ~iscell(loc)
    % HA: vectors/scalar, positive, numeric
    if ~isnumeric(loc) ... % has to be numeric
            || any(isnan(loc)) || any(isinf(loc)) ... % no NaN or Inf
            || ~all(mod(loc,1) == 0) ... % has to be integer
            || ~all(loc > 0) % has to be postive
        throw(CORAerror('CORA:wrongInputInConstructor',...
            ['All entries in the property location have to be '...
            'positive (non-NaN, non-Inf) numeric vectors/scalars.']));
    end
else
    % pHA: each cell-array has to be vectors/scalar, positive, numeric
    for i=1:length(loc)
        if ~isnumeric(loc{i}) ... % has to be numeric
            || any(isnan(loc{i})) || any(isinf(loc{i})) ... % no NaN or Inf
            || ~all(mod(loc{i},1) == 0) ... % has to be integer
            || ~all(loc{i} > 0) % has to be postive
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['All entries in the property location have to be '...
                'positive (non-NaN, non-Inf) numeric vectors/scalars,\n'...
                '  and must not be larger than the number of sets.']));
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
