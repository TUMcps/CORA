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

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   25-April-2023 (MW, disallow R^0)
% Last revision: 10-January-2024 (MW, reformat)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)
    % dimension of space
    dimension = 0;
end

methods

    function obj = fullspace(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1  && isa(varargin{1},'fullspace')
            % copy constructor
            obj = varargin{1};
            return;
        end

        % 2. parse input arguments
        if nargin > 1
            throw(CORAerror('CORA:tooManyInputArgs',1));
        end
        n = varargin{1};

        % 3. check correctness of input arguments
        inputArgsCheck({{n,'att','numeric',{'scalar','nonnegative','integer'}}});
        
        % 4. assign properties
        obj.dimension = n;
    end
end

methods (Static = true)
    fs = generateRandom(varargin) % generate random full space
    fs = enclosePoints(points,varargin) % enclose point cloud with full space
    fs = Inf(n) % instantiates a fullspace fullspace
end

end

% ------------------------------ END OF CODE ------------------------------
