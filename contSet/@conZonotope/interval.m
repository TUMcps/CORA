function I = interval(cZ)
% interval - over-approximate a constrained zonotope object with an
%    axis-aligned interval (bounding box)
%
% Syntax:  
%    I = interval(cZ)
%
% Inputs:
%    cZ - conZonotope object
%
% Outputs:
%    I - interval object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1]; b = 2;
%    cZ = conZonotope(Z,A,b);
%    I = interval(cZ);
%
%    figure; hold on;
%    plot(cZono,[1,2],'FaceColor','r');
%    plot(I,[1,2],'g');
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

if isempty(cZ.A)       % no constraints -> call zonotope method
    
    I = interval(zonotope(cZ.Z));
    
else                    % constraints 
    
    n = size(cZ.Z,1);
    I = interval(zeros(n,1));

    % remove the trivial constraint 0*beta = 0
    cZ = deleteZeros(cZ);

    % loop over all dimensions
    for i = 1:n
        temp = zeros(n,1);
        temp(i) = 1;

        % calculate exact bounds by solving a linear program
        lb = supportFunc(cZ,temp,'lower');
        ub = supportFunc(cZ,temp,'upper');

        I(i) = interval(lb,ub);
    end
end

%------------- END OF CODE --------------