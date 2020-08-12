function check_metasim(options, obj, cDsimRand)
% check_errorOrder - checks if
% options.points|fracVert|fracInpVert|inpChanges (simulateRandom)
% options.points|vertSamp|strechFac (simulateRRT)
%  1) exist
%  2) take an allowed value
%
% Syntax:
%    check_metasim(options, obj)
%
% Inputs:
%    options   - options for object
%    obj       - object of system for case differentiation
%    cDsimRand - if called from simulateRandom
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

% Author:       Mark Wetzlinger
% Written:      08-May-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin < 3
    cDsimRand = false;
end

if cDsimRand
    option = {'points','fracVert','fracInpVert','inpChanges'};
else
    % rrt
    option = {'points','vertSamp','strechFac'};
end
strct = 'options';

% check existence
for i=1:length(option)
    if ~isfield(options,option{i})
        error(printOptionMissing(obj,option{i},strct));
    end
end

% check values
if options.points < 1 || mod(options.points,1) ~= 0
    error(printOptionOutOfRange(obj,'points',strct));
end

if cDsimRand
	if options.fracVert < 0 || options.fracVert > 1
        error(printOptionOutOfRange(obj,'fracVert',strct));
    elseif options.fracInpVert < 0 || options.fracInpVert > 1
        error(printOptionOutOfRange(obj,'fracInpVert',strct));
    elseif options.inpChanges < 0 || mod(options.inpChanges,1) ~= 0
        error(printOptionOutOfRange(obj,'inpChanges',strct));
    end
else
    if ~(options.vertSamp == 0 || options.vertSamp == 1)
        error(printOptionOutOfRange(obj,'vertSamp',strct));
    elseif options.strechFac < 1
        error(printOptionOutOfRange(obj,'strechFac',strct));
    end
end



end

%------------- END OF CODE --------------
