function res = interval(obj)
% interval - over-approximate a constrained zonotope object with an
%            axis-aligned interval (bounding box)
%
% Syntax:  
%    res = interval(obj)
%
% Inputs:
%    obj - c-zonotope object
%
% Outputs:
%    res - interval object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%    int = interval(cZono);
%
%    hold on
%    plotFilled(cZono,[1,2],'r','EdgeColor','none')
%    plot(int,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      13-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if isempty(obj.A)       % no constraints -> call superclass method
    
    res = interval@zonotope(obj);
    
else                    % constraints 
    
    n = size(obj.Z,1);
    res = interval(zeros(n,1));

    % remove the trivial constraint 0*ksi = 0
    obj = removeZeroConstraints(obj);

    % loop over all dimensions
    for i = 1:n
        temp = zeros(n,1);
        temp(i) = 1;

        % calculate exact bounds by solving a linear program
        lb = boundDir(obj,temp,'lower');
        ub = boundDir(obj,temp,'upper');

        res(i) = interval(lb,ub);
    end
end

%------------- END OF CODE --------------