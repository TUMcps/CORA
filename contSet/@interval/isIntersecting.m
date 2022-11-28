function res = isIntersecting(I,S,varargin)
% isIntersecting - determines if an interval intersects a set
%
% Syntax:  
%    res = isIntersecting(I,S)
%    res = isIntersecting(I,S,type)
%
% Inputs:
%    I - interval object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([0;0],[2;2]);
%    I2 = interval([1;1],[3;3]);
%    I3 = interval([-3;-3],[-1;1]);
%
%    isIntersecting(I1,I2)
%    isIntersecting(I1,I3)
%
%    figure; hold on;
%    plot(I1,[1,2],'b');
%    plot(I2,[1,2],'g');
%
%    figure; hold on;
%    plot(I1,[1,2],'b');
%    plot(I3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isIntersecting

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      22-July-2016
% Last update:  14-Sep-2019
%               21-Nov-2019 (NK, added intersection with other sets)
%               12-Mar-2021 (MW, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------

    % pre-processing
    [resFound,vars] = pre_isIntersecting('interval',I,S,varargin{:});

    % check premature exit
    if resFound
        % if result has been found, it is stored in the first entry of var
        res = vars{1}; return
    else
        % assign values
        I = vars{1}; S = vars{2}; type = vars{3};
    end
    
    
    % interval and interval intersection
    if isa(S,'interval')
        
        res = true;
        
        % get object properties
        sup1 = I.sup;
        inf1 = I.inf;
        sup2 = S.sup;
        inf2 = S.inf;
        
        % loop over all dimensions
        for i = 1:length(I)
           if ~isIntersecting1D(inf1(i),sup1(i),inf2(i),sup2(i))
              res = false;
              return
           end
        end
        
    elseif isa(S,'halfspace') || isa(S,'conHyperplane') || ...
           isa(S,'mptPolytope') || isa(S,'ellipsoid')
        
        res = isIntersecting(S,I,type);
        
    else
        
        % exact or over-approximative algorithm
        if strcmp(type,'exact')           
            res = isIntersecting(S,I,type);
        else
            res = isIntersecting(mptPolytope(I),S,type);
        end
    end
    
end


% Auxiliary Functions -----------------------------------------------------

function res = isIntersecting1D(inf1,sup1,inf2,sup2)
% check if two one-dimensional intervals intersect
    res = false;

    if inf1 <= inf2
        if inf2 <= sup1
            res = true;
        end
        
    else % inf2 < inf1
        if inf1 <= sup2
            res = true;
        end
        
    end

end

%------------- END OF CODE --------------