classdef transition
% transition - constructor of the class transition
%
% Syntax:  
%    obj = transition(guard,reset,target)
%
% Inputs:
%    guard - guard set specified as contSet object
%    reset - reset function (only linear map!): Ax+b, with struct:
%            reset.A, reset.b
%    target - target: int (id of target location)
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
end

methods
    
    % class constructor
    function obj = transition(varargin)

        % parse input arguments
        if nargin == 3
            
            % assign object properties
            obj.guard = varargin{1};
            obj.reset = varargin{2};
            obj.target = varargin{3};
            
            % check reset function
            list = fields(obj.reset);
            
            if isfield(obj.reset,'f')
                tmp = list(cellfun(@(x) ~strcmp(x,'f'),list));
                if ~isempty(tmp)
                    tmp = cellfun(@(x) [x,','],tmp,'UniformOutput',false);
                    tmp = [tmp{:}]; tmp = tmp(1:end-1);
                    warning(['The following fields of struct "reset"', ...
                                                 ' are redundant: ',tmp]);
                end
            else
                if ~isfield(obj.reset,'A')
                    error('Field "A" of struct "reset" is missing!');
                elseif ~isfield(obj.reset,'b')
                    error('Field "b" of struct "reset" is missing!');
                else
                    tmp = list(cellfun(@(x) ~strcmp(x,'A') && ...
                                                     ~strcmp(x,'b'),list));
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
      
        n = dim(obj.guard); x = sym('x',[n,1]);
        f = obj.reset.f(x);
        
        % first order derivative
        J = jacobian(f,x);
        obj.reset.J = matlabFunction(J,'Vars',{x});
        
        % second order derivative
        obj.reset.Q = cell(n,1);
        
        for i = 1:n
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