function options = check_timeStepLoc(options, obj)
% check_timeStepLoc - checks if options.timeStepLoc
%  1) exists
%  2) takes an allowed value
%
% Syntax:
%    options = check_timeStepLoc(options, obj)
%
% Inputs:
%    options - options for object
%    obj     - system object
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      08-May-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    option = 'timeStep';
    strct = 'options';

    % check if timeStepLoc is specified -> warning that it is redundant
    if isfield(options,'timeStepLoc')
        warning('options.timeStepLoc is redundant!');
    end

    % get number of locations
    locations = obj.location;
    numLoc = length(locations);

    % assign internal value options.Uloc which stores the input set for each
    % location
    if ~isfield(options,option)
    % display error message
        error(printOptionMissing(obj,option,strct));      

    elseif isscalar(options.timeStep) 
    % same input set for each location    

        timeStepLoc = cell(numLoc,1);

        for i = 1:numLoc
            timeStepLoc{i} = options.timeStep;
        end

        options = rmfield(options,option);

    elseif iscell(options.timeStep)
    % copy time step for each location

        if all(size(options.timeStep) ~= [numLoc,1]) && ...
           all(size(options.timeStep) ~= [1,numLoc])
            error(printOptionOutOfRange(obj, option, strct));
        else
            timeStepLoc = options.timeStep;
        end

        options = rmfield(options,option);

    else
        error(printOptionOutOfRange(obj, option, strct));
    end

    options.timeStepLoc = timeStepLoc;

end

%------------- END OF CODE --------------