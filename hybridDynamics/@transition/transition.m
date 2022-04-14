classdef transition
% transition - constructor of the class transition
%
% Syntax:  
%    obj = transition(guard,reset,target)
%
% Inputs:
%    guard - guard set specified as contSet object
%    reset - reset function (linear map or nonlinear function handle)
%    target - target: int (id of target location)
%    stateDim - dimension of state space (in case it is not discernable by 
%               other arguments)
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton, location

% Author:       Matthias Althoff
% Written:      02-May-2007 
% Last update:  30-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    guard = [];         % guard set
    reset = [];         % reset function 
    target = [];        % target location
    synchLabel = [];    % synchronization label
    stateDim;           % state space dimension
end

methods
    
    % class constructor
    function obj = transition(varargin)

        % parse input arguments
        if nargin >= 3 && nargin <= 5
            
            % assign object properties
            obj.guard = varargin{1};
            obj.reset = varargin{2};
            obj.target = varargin{3};
            obj.synchLabel = [];
            if nargin >= 4
                obj.synchLabel = varargin{4};
            end
            
            if isempty(obj.guard) && ~isfield(obj.reset,'A') && nargin < 5
                error("Transitions with empty guard! In this case, also specify "...
                    + "the state space dimension as constructor argument!");
            end
            
            if nargin == 5
                obj.stateDim = varargin{5};
            end
            
            % check reset function
            list = fields(obj.reset);
            
            
            resetInput = 0;
            % check if reset function is dependent on inputs
            if isfield(obj.reset,'hasInput')
                resetInput = obj.reset.hasInput;
            end
           
            if isfield(obj.reset,'f')
                
                % for nonlinear reset functions which depend on inputs,
                % there is no way to determine the input dimension from the
                % transition object, therefore it has to be a field of the
                % reset struct
                if resetInput && ~isfield(obj.reset,'inputDim')
                    error('Field "inputDim" of struct "reset" is missing!');
                end
                   
                tmp = list(cellfun(@(x) ~strcmp(x,'f') && ~strcmp(x,'hasInput') ...
                                        && ~strcmp(x,'inputDim') ...
                                        && ~strcmp(x,'stateDim'),list));
                if ~isempty(tmp)
                    tmp = cellfun(@(x) [x,','],tmp,'UniformOutput',false);
                    tmp = [tmp{:}]; tmp = tmp(1:end-1);
                    warning(['The following fields of struct "reset"', ...
                                                 ' are redundant: ',tmp]);
                end
            else
                if ~isfield(obj.reset,'A')
                    error('Field "A" of struct "reset" is missing!');
                elseif ~isfield(obj.reset,'c')
                    error('Field "c" of struct "reset" is missing!');
                % if the resut function depends on inputs, we need a B
                % matrix
                elseif resetInput && ~isfield(obj.reset,'B')
                    error('Field "B" of struct "reset" is missing!');
                else
                    tmp = list(cellfun(@(x) ~strcmp(x,'A') && ~strcmp(x,'hasInput')...
                        && ~strcmp(x,'inputDim') && ~strcmp(x,'c') ...
                        && ~strcmp(x,'B') && ~strcmp(x,'stateDim'),list));
                    if ~isempty(tmp)
                        tmp = cellfun(@(x) [x,','],tmp, ...
                                                    'UniformOutput',false);
                        tmp = [tmp{:}]; tmp = tmp(1:end-1);
                        warning(['The following fields of struct' ... 
                                          ' "reset" are redundant: ',tmp]);
                    end
                end
            end
            
            % pre-compute derivaties for nonlinear resets
            if isfield(obj.reset,'f')
                obj = compDerivatives(obj);
            end
        else
            error('Wrong number of input arguments!');
        end
    end
    
    % auxiliary functions
    function obj = compDerivatives(obj)
    % precompute all derivatives required for the enclosure of the 
    % nonlinear reset function by a Taylor series expansion of order 2   
      
    
        % check if reset function is dependent on inputs
        resetInput = 0;
        if isfield(obj.reset,'hasInput')
            resetInput = obj.reset.hasInput;
        end
        
        % state dimension n
        if isempty(obj.guard)
            n = obj.stateDim;
        else
            n = dim(obj.guard);
        end
        
        % create symbolic variables for reset function
        x = sym('x',[n,1]);
        
        % input dimension m
        m = 0;
        % if reset function depends on inputs, we create the subsitute
        % state x' = [x;u] which is used for the rest of this function
        if resetInput
            m = obj.reset.inputDim;
            u = sym('u',[m,1]);
            
            % substitute state x for extended state x' = [x;u]
            x = [x;u];
            % alter reset function to accept singular input x'
            obj.reset.f = @(x) obj.reset.f(x(1:n),x(n+1:end));
            % alter state dimension according to dimension of x'
            n = n+m;
        end
        
        % evaluta f symbolically
        f = obj.reset.f(x);
        
        % first order derivative
        J = jacobian(f,x);
        obj.reset.J = matlabFunction(J,'Vars',{x});
        
        % second order derivative
        obj.reset.Q = cell(n,1);
        
        for i = 1:n-m
            obj.reset.Q{i,1} = matlabFunction(hessian(f(i),x),'Vars',{x});
        end
        
        % third order derivative
        obj.reset.T = cell(size(J));
        
        for i = 1:size(J,1)
            for j = 1:size(J,2)
                temp = hessian(J(i,j),x);
                if any(any(temp ~= 0))
                    obj.reset.T{i,j} = matlabFunction(temp,'Vars',{x});
                end
            end
        end
    end

end
end

%------------- END OF CODE --------------