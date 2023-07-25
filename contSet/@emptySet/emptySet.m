classdef emptySet < contSet
% emptySet - object constructor for empty sets
%
% Description:
%    This class represents empty sets.
%
% Syntax:
%    obj = emptySet()
%    obj = emptySet(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    obj - generated emptySet object
%
% Example: 
%    n = 2;
%    O = emptySet(n);
%    plot(O);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Mark Wetzlinger
% Written:      22-March-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    % dimension of space
    dimension = 0;
end

methods

    function obj = emptySet(varargin)
        
        if nargin == 0
            % nothing to see here...

        elseif nargin == 1 
            if isa(varargin{1},'emptySet')
                % copy constructor
                obj = varargin{1};
                return;
            else
                % instantiate n-dimensional space
                n = varargin{1};
                inputArgsCheck({{n,'att','numeric',{'scalar','nonnegative','integer'}}});
                
                % set property
                obj.dimension = n;
            end

        else
            throw(CORAerror('CORA:tooManyInputArgs',1));
        end

    end
end

methods (Static = true)
    O = generateRandom(varargin) % generate random empty set
    O = enclosePoints(points,varargin) % enclose point cloud with empty
end

end

%------------- END OF CODE --------------