function options = check_tStart_tFinal(obj, options, checkName)
% check_tStart_tFinal - checks if options.tStart/tFinal
%  1) exist
%  2) take an allowed value
%
% Syntax:
%    check_tStart_tFinal(obj, options, checkName)
%
% Inputs:
%    obj       - system object
%    options   - options for object
%    checkName - options check (reach or simulate)
%
% Outputs:
%    options - options for object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Mark Wetzlinger
% Written:      04-Mar-2019
% Last update:  08-Aug-2019
%               03-May-2020 (rewriting of error msgs using class(obj))
% Last revision:---

%------------- BEGIN CODE --------------

strct = 'params';

option1 = 'tStart';
defValue = getDefaultOption(option1);
% check tStart ... optional, default at 0
if ~isfield(options,option1)
    options.tStart = defValue;
elseif options.tStart == defValue
    % not necessary since default value
    printDefaultValue(obj,option1,defValue);
elseif options.tStart < 0
    error(printOptionOutOfRange(obj,option1,strct));
end

option2 = 'tFinal';
% check tFinal ... mandatory
if ~isfield(options,option2)
    error(printOptionMissing(obj,option2,strct));
elseif options.tFinal < options.tStart
    error(printOptionOutOfRange(obj,option2,strct));
elseif isa(obj,'linearSysDT') || isa(obj,'nonlinearSysDT') % discrete-time systems
    t = options.tStart:obj.dt:options.tFinal;
    if t(end) ~= options.tFinal
        error(printOptionOutOfRange(obj,option2,strct));
    end
end

% check if final time unnecessary large
% only for linear systems, not for simulation
if ~strcmp(checkName, 'checkOptionsSimulate') && isa(obj,'linearSys') && any(obj.A(:))
    try
        if ~issparse(obj.A)
            maxEig = max(real(eig(obj.A)));
        else
            maxEig = real(eigs(obj.A,1,'smallestabs'));
        end
        if maxEig < 0
            threshold_tFinal = 0.001;
            tSettled = log(threshold_tFinal)/maxEig;
            if options.tFinal > tSettled
                options.tFinal = tSettled;
            end
        end
    catch
        % keep options.tFinal as it is, if A singular
    end
end


end

%------------- END OF CODE --------------
