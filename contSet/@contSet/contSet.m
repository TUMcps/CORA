classdef contSet
% contSet - abstract superclass for continuous sets
%
% Syntax:  
%    S = contSet()
%    S = contSet(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    obj - generated contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       02-May-2007 
% Last update:   04-May-2020 (MW, transition to classdef)
%                01-June-2022 (MW, add CORAerror)
%                22-March-2023 (MW, remove deprecated property dimension)
% Last revision: 16-June-2023 (MW, restructure using standardized workflow)

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    % no properties
end

methods

    function obj = contSet(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'contSet')
            obj = varargin{1}; return
        end
        
        % 2. parse input arguments
        if nargin > 0
            throw(CORAerror('CORA:tooManyInputArgs',1));
        end

        % 3. assign properties
        % (none)
    end
end

methods (Static = true)
    S = generateRandom(varargin) % generates a random contSet
end

end

%------------- END OF CODE --------------