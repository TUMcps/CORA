function obj = plus(obj,S)
% plus - Overloaded '+' operator for the Minkowski addition of a set or a
%        vector with a reachSet object
%
% Syntax:  
%    obj = plus(obj,S)
%
% Inputs:
%    obj - reachSet object 
%    S - contSet object or numerical vector
%
% Outputs:
%    obj - resulting tranformed reachset object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Niklas Kochdumper
% Written:      04-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % get reachSet object
    if ~isa(obj,'reachSet')
       temp = obj;
       obj = S;
       S = temp;
    end

    % compute Minkowski sum
    for i = 1:size(obj,1)
       obj(i).timeInterval.set = cellfun(@(x) S + x, ...
                            obj(i).timeInterval.set,'UniformOutput',false);
       obj(i).timePoint.set = cellfun(@(x) S + x, ...
                               obj(i).timePoint.set,'UniformOutput',false);             
    end

%------------- END OF CODE --------------