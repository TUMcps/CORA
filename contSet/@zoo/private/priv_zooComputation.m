function res = priv_zooComputation(fhandle, obj1, obj2)
% priv_zooComputation - executive function
%
% Syntax:
%    res = priv_zooComputation(fhandle, obj1, obj2)
%
% Inputs:
%    fhandle - function handle
%    obj1 - a zoo object
%    obj2 - a zoo object
%
% Outputs:
%    res - a zoo object
%
% Example: 
%
% Other m-files required: taylm, interval
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, interval

% Authors:       Dmitry Grebenyuk
% Written:       11-November-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    if nargin == 2

        res = arrayfun(@(a) aux_s_zooCmp1(fhandle, a), obj1, 'UniformOutput', 0);

    elseif nargin == 3

        res = arrayfun(@(a) aux_s_zooCmp2(fhandle, a, obj2), obj1, 'UniformOutput', 0);

    end
    A = cat(1, res{:});
    res = reshape(A, size(res));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_s_zooCmp1(fhandle, obj)
   
    res = obj;
    
    for i = 1:length(obj.method)
       res.objects{i} = fhandle(obj.objects{i}); 
    end
end

function res = aux_s_zooCmp2(fhandle, obj1, obj2)
   
    res = obj1;
    
    for i = 1:length(obj1.method)
       res.objects{i} = fhandle(obj1.objects{i},obj2); 
    end

end

% ------------------------------ END OF CODE ------------------------------
