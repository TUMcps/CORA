classdef hybridAutomaton
% hybridAutomaton - constructor for class hybridAutomaton
%
% Syntax:  
%    obj = hybridAutomaton(loc)
%
% Inputs:
%    loc - cell-array storing the location objects
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton, location, transition

% Author: Matthias Althoff
% Written: 03-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    location = [];       % cell-array of location objects
end

methods
    
    % class constructor
    function obj = hybridAutomaton(varargin)

        % parse input arguments
        if nargin == 1
            obj.location = varargin{1};
        else
            error('Wrong number of input arguments!');
        end
    end
end
end

%------------- END OF CODE --------------