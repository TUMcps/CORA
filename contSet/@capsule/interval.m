function res = interval(obj)
% interval - Over-approximate a capsule by an interval
%
% Syntax:  
%    res = interval(obj)
%
% Inputs:
%    obj - capsule object
%
% Outputs:
%    res - interval object
%
% Example: 
%    C = capsule([0;0],[-2;2],2);
%
%    int = interval(C);
%
%    figure
%    hold on
%    plot(C,[1,2],'k');
%    plot(int,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Author:       Niklas Kochdumper
% Written:      20-Nov-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % initialization
    n = length(obj.c);
    sup = zeros(n,1);
    infi = zeros(n,1);
    E = eye(n);
    
    % loop over all dimensions
    for i = 1:n
        sup(i) = supportFunc(obj,E(:,i),'upper');
        infi(i) = supportFunc(obj,E(:,i),'lower');        
    end
    
    % construct resulting interval
    res = interval(infi,sup);

%------------- END OF CODE --------------