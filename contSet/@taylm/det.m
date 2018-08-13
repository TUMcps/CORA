function res = det(obj)
% det - Calculate the determinant of a taylor model matrix
%
% Syntax:  
%    res = det(obj)
%
% Inputs:
%    obj - a taylm object
% Outputs:
%    res - a taylm object
%
% Example: 
%    t = taylm(interval(1,2),10,'x');
%    T = [sin(t) 2*t; t exp(t)];
%    D = det(T)
%
% Other m-files required: taylm
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, traace

% Author:       Niklas Kochdumper
% Written:      12-June-2018
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

    % get matrix size
    [~,n] = size(obj);
    
    % trivial case
    if n == 1
        res = obj;
        
    % recursive function call
    else
        res = 0;
        for j = 1:n
            A1 = obj(2:n, [1:j-1, j+1:n]);
            res = res + (-1)^(j+1)*obj(1,j)*det(A1);
        end
    end
    
%------------- END OF CODE --------------