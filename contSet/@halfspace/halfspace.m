classdef halfspace < contSet
% halfspace - object constructor for halfspaces 
%
% Description:
%    This class represents halfspace objects defined as
%    {x | c'*x <= d}.
%
% Syntax:
%    obj = halfspace(hs)
%    obj = halfspace(c,d)
%
% Inputs:
%    c - normal vector of the halfspace (n-by-1)
%    d - distance to the origin (scalar)
%
% Outputs:
%    obj - generated halfspace object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: example_halfspace.m

% Authors:       Matthias Althoff
% Written:       06-June-2011
% Last update:   14-August-2019
%                02-May-2020 (add property validation)
%                19-March-2021 (MW, error messages)
%                14-December-2022 (TL, property check in inputArgsCheck)
%                24-July-2023 (MW, disable 0*x <= d)
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------


properties (SetAccess = private, GetAccess = public)
    c = [];     % normal vector
    d = 0;      % offset
end
    
methods
    %class constructor
    function obj = halfspace(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'halfspace')
            obj = varargin{1}; return
        end
        
        % 1. parse input arguments: varargin -> vars
        [c,d] = aux_parseInputArgs(varargin{:});

        % 2. check correctness of input arguments
        aux_checkInputArgs(c,d,nargin);

        % 3. assign properties (reshape c to a column vector)
        obj.c = reshape(c,[],1);
        obj.d = d;
    end
end

methods (Static = true)
    hs = generateRandom(varargin) % generates random halfspace
    hs = empty(n) % instantiates an empty halfspace
    hs = Inf(n) % instantiates a fullspace halfspace
end

end


% Auxiliary functions -----------------------------------------------------

function [c,d] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 2
        throw(CORAerror('CORA:tooManyInputArgs',2));
    end

    % set default values
    [c,d] = setDefaultValues({[],[]},varargin);

end

function aux_checkInputArgs(c,d,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        inputArgsCheck({ ...
            {c, 'att', 'numeric', {'finite', 'vector'}}; ...
            {d, 'att', 'numeric', {'finite', 'scalar'}}; ...
        })

        % ensure that halfspace actually constrains some value
        % ...if this functionality is desired, use fullspace class instead
%         if ~any(c)
%             throw(CORAerror('CORA:wrongInputInConstructor',...
%                 'At least one entry in the normal vector has to be non-zero.'));
%         end
        
    end

end

% ------------------------------ END OF CODE ------------------------------
