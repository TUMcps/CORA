function obj = mtimes(M,obj)
% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%          with a reachSet object
%
% Syntax:  
%    obj = mtimes(M,obj)
%
% Inputs:
%    M - numerical matrix
%    obj - reachSet object 
%
% Outputs:
%    obj - resulting tranformed reachset object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Niklas Kochdumper
% Written:      04-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    for i = 1:size(obj,1)
       obj(i).timeInterval.set = cellfun(@(x) M*x, ...
                            obj(i).timeInterval.set,'UniformOutput',false);
       obj(i).timePoint.set = cellfun(@(x) M*x, ...
                               obj(i).timePoint.set,'UniformOutput',false);             
    end

%------------- END OF CODE --------------