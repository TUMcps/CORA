function I = interval(z)
% interval - calculate the bounding interval of a zoo-object
%
% Syntax:
%    I = interval(z)
%
% Inputs:
%    z - zoo object
%
% Outputs:
%    I - interval overapproximating the zoo (class interval)
%
% Other m-files required: interval
% Subfunctions: none
% MAT-files required: none
%
% See also: zoo

% Authors:       Niklas Kochdumper
% Written:       10-April-2018
% Last update:   --- 
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    I = arrayfun(@(a) aux_zoo2int(a), z, 'UniformOutput', false);
    A = [I{:}];
    I = reshape(A, size(I)); 
    
end


% Auxiliary functions -----------------------------------------------------

function res = aux_zoo2int(z)

    res = interval(-Inf,Inf);
    for i = 1:length(z.method)
        res = and_(res,interval(z.objects{i}),'exact'); 
    end
    
end
    
% ------------------------------ END OF CODE ------------------------------
