classdef fullspace < contSet
% fullspace - object constructor for n-dimensional spaces; setting n=0 is
%    only permitted for some functions, because the 0-dimensional 0 vector
%    cannot be represented in MATLAB (different from the 1-dimensional
%    '0'). Still, the results of the remaining set operations are
%    described in the respective function headers.
%
% Description:
%    This class represents objects defined as {x \in \R{n}}.
%
% Syntax:
%    obj = fullspace()
%    obj = fullspace(n)
%
% Inputs:
%    n - dimension
%
% Outputs:
%    obj - generated fullspace object
%
% Example: 
%    n = 2;
%    fs = fullspace(n);
%    plot(fs);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: emptySet

% Author:       Mark Wetzlinger
% Written:      22-March-2023
% Last update:  25-April-2023 (MW, disallow R^0)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    % dimension of space
    dimension = 0;
end

methods

    function obj = fullspace(varargin)
        % number of input arguments
        if nargin > 1
            throw(CORAerror('CORA:tooManyInputArgs',1));
        end
        
        if nargin == 0
            % nothing to see here...

        elseif nargin == 1 
            if isa(varargin{1},'fullspace')
                % copy constructor
                obj = varargin{1};
                return;
            else
                % instantiate n-dimensional space
                n = varargin{1};
                inputArgsCheck({{n,'att','numeric',...
                    {'scalar','nonnegative','integer'}}});
                
                % set property
                obj.dimension = n;
            end
        end

    end
end

methods (Static = true)
    fs = generateRandom(varargin) % generate random full space
    fs = enclosePoints(points,varargin) % enclose point cloud with full space
end

end

%------------- END OF CODE --------------