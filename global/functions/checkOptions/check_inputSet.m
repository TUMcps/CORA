function options = check_inputSet(obj, options)
% check_inputSet - checks if options.U|u
%  1) exist
%  2) take an allowed value
%
% Syntax:
%    check_inputSet(obj, options)
%
% Inputs:
%    obj       - linear system
%    options   - options for object
%
% Outputs:
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: set_inputSet

% Author:       Mark Wetzlinger
% Written:      04-Mar-2019
% Last update:  18-Dec-2019
%               03-May-2020 (init empty U, error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'params';
% input sizes must match
if ~isfield(options,'U')
    if ~(isa(obj,'hybridAutomaton') || isa(obj,'parallelHybridAutomaton'))
        options.U = zonotope(zeros(max(obj.nrOfInputs,1),1));
    end
elseif isa(obj,'linearSys')
    if isscalar(obj.B)
        if length(center(options.U)) ~= obj.dim
            % if obj.B is scalar, U has to have same dim as obj.A
            error(printOptionOutOfRange(obj, 'U', strct));
        end
    elseif size(obj.B,2) ~= length(center(options.U))
        % U has to be compatible with #cols of B (if B ismatrix)
        error(printOptionOutOfRange(obj, 'U', strct));
    end
end
    

if isfield(options,'u')
    % u is given
    option = 'u';
    if isfield(options,'linAlg') && strcmp(options.linAlg,'adap')
        % u cannot be given when adaptive method used
        error(printOptionSpecificError(obj,option,...
            'No u when options.linAlg = adap.'));
    elseif size(center(options.U),1) ~= size(options.u,1)
        % u has to have same dim as U
        error(printOptionSpecificError(obj,option,...
            'U and u must have the same dimension.'));
    elseif any(center(options.U))
        % U has to be centered at origin
        error(printOptionOutOfRange(obj, 'U', strct));
    elseif isfield(options,'timeStep') && size(options.u,2) ~= ...
           ceil((options.tFinal - options.tStart)/options.timeStep)
        % u has to equal the number of time steps
        error(printOptionOutOfRange(obj, 'u', strct));
    elseif (isa(obj,'linearSysDT') || isa(obj,'nonlinearSysDT')) && ...
            size(options.u,2) ~= ceil((options.tFinal - options.tStart)/obj.dt)+1
        % u has to equal the number of time steps
        error(printOptionOutOfRange(obj, 'u', strct));
    end
end

end

%------------- END OF CODE --------------
