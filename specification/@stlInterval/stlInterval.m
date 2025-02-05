classdef stlInterval < contSet
% stlInterval - open/closed time intervals for signal temporal logic specifications
%
% Time intervals are always one dimensional and are subsets of the non-negative reals.
% The boundaries of the interval can be open or closed. The interval can be empty.
%
% Syntax:
%    int = stlInterval()
%    int = stlInterval(x)
%    int = stlInterval(I)
%    int = stlInterval(SI)
%    int = stlInterval(lb,ub)
%    int = stlInterval(lb,ub,closed)
%    int = stlInterval(lb,ub,lc,rc)
%
% Inputs:
%    x - numeric value
%    I - interval object
%    SI - stlInterval object
%    lb - lower bound
%    ub - upper bound
%    closed - boolean for both left and right closed (default: true)
%    lc - is left closed (default: true)
%    rc - is right closed (default: true)
%
% Outputs:
%    int - generated stlInterval object
%
% Example:
%
%    int0 = stlInterval();               % empty interval
%    int1 = stlInterval(1);              % singular interval [1,1]
%    int2 = stlInterval(0,1);            % closed interval [0,1]
%    int3 = stlInterval(0,1,false);      % open interval (0,1)
%    int4 = stlInterval(0,1,true,false); % half open interval [0,1)
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: 15-October-2024 (MW, restructure)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    lower {mustBeScalarOrEmpty,mustBeNumeric}
    leftClosed logical = true
    upper {mustBeScalarOrEmpty,mustBeNumeric}
    rightClosed logical = true
end

methods

    % class constructor
    % TODO: maybe add possibility to create from interval and convert all intervals in STL to this class
    function obj = stlInterval(varargin)

        % 0. check number of input arguments
        assertNarginConstructor(0:4,nargin);

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'stlInterval')
            obj = varargin{1};
            return
        end

        % 2. parse input arguments: varargin -> vars
        [lb,ub,lc,rc] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(lb,ub,lc,rc,nargin);

        % 4. compute properties
        [lb,ub,lc,rc] = aux_computeProperties(lb,ub,lc,rc);

        % 5. assign properties
        obj.lower = lb;
        obj.upper = ub;
        obj.leftClosed = lc;
        obj.rightClosed = rc;

        % 6. set precedence (via contSet)
        obj.precedence = 13;
    end

    % methods in separate files
    n = dim(I)
    [infi,isMin] = infimum(obj)
    res = isequal(obj,other)
    int = leftClosure(obj)
    int = minkDiff(obj,sub,type)
    int = minus(minuend,subtrahend)
    int = plus(summand1,summand2)
    int = rightClosure(obj)
    [sup,isMax] = supremum(obj)
    int = toLeft(obj)
    int = toRight(obj)

    % overload equality and inequality operators using isequal
    function res = eq(obj,other)
        res = isequal(obj,other);
    end
    function res = ne(obj,other)
        res = ~isequal(obj,other);
    end

    function str = toStr(obj)
        if isemptyobject(obj)
            str = '∅';
            return;
        end

        if obj.leftClosed
            lterm = '[';
        else
            lterm = '(';
        end
        if obj.rightClosed
            rterm = ']';
        else
            rterm = ')';
        end
        str = sprintf('%s%d, %d%s',lterm,obj.lower,obj.upper,rterm);
    end

    % override horzcat and vertcat
    function out = horzcat(varargin)
        out = builtin('horzcat',varargin{:});
    end
    function out = vertcat(varargin)
        out = builtin('vertcat',varargin{:});
    end
    function out = cat(varargin)
        out = builtin('cat',varargin{:});
    end
    function out = subsasgn(varargin)
        out = builtin('subsasgn',varargin{:});
    end
end

methods (Static = true)
    I = empty(n)
end
end


% Auxiliary functions -----------------------------------------------------

function [lb,ub,lc,rc] = aux_parseInputArgs(varargin)
% parse inputs
    
    lb = []; ub = []; lc = true; rc = true;

    [lb,ub,lc,rc] = setDefaultValues({lb,ub,lc,rc},varargin);

    if nargin == 1
        ub = lb;
    elseif nargin == 3
        rc = lc;
    end
    
end

function aux_checkInputArgs(lb,ub,lc,rc,n_in)
% check inputs

if CHECKS_ENABLED && n_in > 0
    if ~isempty(lb) && lb < 0
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Lower bound must be non-negative.'));
    end
end

end

function [lb,ub,lc,rc] = aux_computeProperties(lb,ub,lc,rc)
% compute properties

    if ~isempty(lb) && ~isempty(ub)
        if abs(lb) == inf
            lc = false;
        end
        if abs(ub) == inf
            rc = false;
        end
        if lb > ub || (lb == ub && (~lc || ~rc))
            % make empty interval
            [lb,ub,lc,rc] = aux_parseInputArgs();
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
