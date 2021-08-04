function res = in(obj1,obj2,varargin)
% in - determines if obj2 is contained in obj1
%
% Syntax:  
%    res = in(obj1,obj2)
%    res = in(obj1,p,tol)
%
% Inputs:
%    obj1 - mptPolytope object
%    obj2 - conSet object
%    p - single point
%    tol - numerical tolerance for point in set containment
%
% Outputs:
%    res - 1/0 if obj2 is contained in obj1, or not
%
% Example: 
%    poly1 = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    poly2 = mptPolytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]);
%    poly3 = poly2 + [2;0];
%
%    in(poly1,poly2)
%    in(poly1,poly3)
%
%    figure
%    hold on
%    plot(poly1,[1,2],'b');
%    plot(poly2,[1,2],'g');
%
%    figure
%    hold on
%    plot(poly1,[1,2],'b');
%    plot(poly3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/in

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  26-July-2021 (VG: extended to multiple points)
% Last revision:---

%------------- BEGIN CODE --------------

    res = 0;
    
    % get object properties
    A = obj1.P.A;
    b = obj1.P.b;
    
    % parse input arguments
    tol = 0;

    if nargin >= 3 && ~isempty(varargin{1})
       tol = varargin{1}; 
    end
    
    % point in polytope containment
    if isnumeric(obj2)
        
        if all(all(A*obj2 - b <= tol))
            res = 1;
        end

    % other set in polytope containment
    else

        % loop over all halfspaces
        for i = 1:size(A,1)
            b_ = supportFunc(obj2,A(i,:)','upper');
            if b_ > b(i)+tol
               return 
            end
        end
        
        res = 1;

    end
end

%------------- END OF CODE --------------
