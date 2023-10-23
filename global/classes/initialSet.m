classdef initialSet
% initialSet - class that stores the initial set of a reachSet object
%
% Syntax:
%    R0 = initialSet(set)
%
% Inputs:
%    set - initial set
%
% Outputs:
%    R0 - generated initialSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: reach

% Authors:       Tobias Ladner
% Written:       01-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
    properties
        set
    end
    
    methods
        % constructor
        function obj = initialSet(set)
            % parse input
            if nargin < 1
                throw(CORAerror('CORA:notEnoughInputArgs', 1))
            elseif nargin > 1
                throw(CORAerror('CORA:tooManyInputArgs', 1))
            end
            inputArgsCheck({{set, 'att', 'contSet'}})
            obj.set = set;
        end
        
        function han = plot(R0, varargin)
            dims = setDefaultValues({[1,2]},varargin);
            NVpairs = readPlotOptions(varargin(2:end),'initialSet');
            han = plot(R0.set, dims, NVpairs{:});

            if nargout == 0
                clear han;
            end
        end
        
        function han = plotOverTime(R0, varargin)
            dims = setDefaultValues({1},varargin);
            inputArgsCheck({{R0,'att','initialSet'};
                    {dims,'att',{'numeric'},{'nonempty','scalar','integer','positive'}}});
            NVpairs = readPlotOptions(varargin(2:end),'initialSet');
    
            set = project(R0.set, dims);
            set = [0; 1] * set;
            han = plot(set, [1, 2], NVpairs{:});

            if nargout == 0
                clear han;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
