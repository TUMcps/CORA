function res = isIntersecting(obj1,obj2,option)
% isIntersecting - determines if elements of a zonotope obj2
% are in another zonotope obj1 -> approximative intersection calculated
%
% Syntax:  
%    res = isIntersecting(obj1,obj2,option)
%
% Inputs:
%    obj1 - 1st zonotope object
%    obj2 - 2nd zonotope object
%    option - none or 'approx'
%
% Outputs:
%    res  - 1/0 if zonotope intersects, or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:        Matthias Althoff
% Written:       07-May-2007 
% Last update:   14-Sep-2019 (MW, rename in -> isIntersecting)
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin == 3
    if ~strcmp(option,'approx')
        error("Wrong flag in isIntersecting");
    end
    approx = true;
else
    warning("Overapproximative intersection computed");
    approx = false;
end


if approx
    
    %simple test: Is the center of obj2 in obj1?
    c=center(obj2);
    inequality=(obj1.halfspace.H*c<=obj1.halfspace.K);

    if all(inequality)
        res=1;
    end

    %use idea of the previous interval test
    %make a polytope of obj1
    P=polytope(obj1.halfspace.H,obj1.halfspace.K);

    %test 1: Is the center of obj2 in obj1?
    c=center(obj2);
    bool=isinside(P,c);

    if bool==1
        res=1;
    else
        %test 2: Is any vertex of obj2 in obj1?
        %overapproximate zonotope to order 1 to prevent an explosion of the
        %number of vertices
        obj2=reduce(obj2,'girard',1);
        %get vertices of the zootope
        V=vertices(obj2);
        %check if vertex in polytope
        bool=0;
        iVertex=1;
        while (bool==0) & (iVertex<=length(V(:,1)))
            bool=isinside(P,V(iVertex,:));
            iVertex=iVertex+1;
        end
        if bool==1
            res=1;
        else
            %test 3: check if the intersection of obj1 and obj2 is empty
            intersection=P&polytope(obj2);
            if ~isempty(intersection)
                res=1;
            else
                res=0;
            end
        end
    end
    
else
    error("This flag is not valid/yet implemented");
end


%------------- END OF CODE --------------