classdef (Abstract) logicSignal
% logicSignal - interface for signals of logic values (truth values) over time
%
% Syntax:
%    logicSignal()
%
% Inputs:
%    none (class is abstract)
%
% Outputs:
%    none (class is abstract)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: pointSegmentSignal, kleeneSignal, fourValuedSignal, finiteSignal

% Authors:       Florian Lercher
% Written:       09-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

methods (Abstract)
    % basic operations
    val = at(obj,time) % get signal value at time
    sig = set(obj,interval,value) % set signal value in interval
    
    % logical operations
    sig = not(obj) % logical not
    sig = until(lhs,interval,rhs) % logical until
    
    plot(obj) % plot the signal
end

methods
    % logical and
    function sig = and(obj,other)
        sig = obj.and_(obj,other);
    end

    % logical or
    function sig = or(obj,other)
        sig = obj.or_(obj,other);
    end
end

methods (Abstract, Static)
    % methods to combine multiple signals point-wise
    sig = and_(varargin) % logical and
    sig = or_(varargin) % logical or
    sigs = combine(op,varargin) % arbitrary point-wise operation op
end
end

% ------------------------------ END OF CODE ------------------------------
