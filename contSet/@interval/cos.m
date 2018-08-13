function res = cos(intVal)
% cos - Overloaded 'cos()' operator for intervals
%
% Syntax:  
%    res = cos(intVal)
%
% Inputs:
%    intVal - interval object
%
% Outputs:
%    res - interval object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Matthias Althoff, Dmitry Grebenyuk
% Written:      25-June-2015
% Last update:  06-January-2016 (DG)
%               05-February-2016 (MA)
% Last update:  22-February-2016 the matrix case is rewritten  (Dmitry Grebenyuk)
% Last revision:---

% inf is x infimum, sup is x supremum

% [-1, 1]                       if (sup - inf) >= 2*pi,
% [-1, 1]                       if (sup - inf) < 2*pi and inf <= pi and sup <= pi and sup < inf,
% [cos(sup), cos(inf)]          if (sup - inf) < 2*pi and inf <= pi and sup <= pi and sup >= inf,
% [-1, max(cos(inf),cos(sup))]  if (sup - inf) < 2*pi and inf <= pi and sup > pi,
% [-1, 1]                       if (sup - inf) < 2*pi and inf > pi and sup > pi and sup < inf,
% [min(cos(inf),cos(sup)), 1]   if (sup - inf) < 2*pi and inf > pi and sup <= pi,
% [cos(inf), cos(sup)]          if (sup - inf) < 2*pi and inf > pi and sup > pi and sup >= inf.

%------------- BEGIN CODE --------------

% scalar case
if isnumeric(intVal)
    
        res = interval();

        %sup - inf >= 2pi
        if intVal.sup - intVal.inf >= 2*pi
            res.inf = -1;
            res.sup = 1;
        else
            %remove multiples of 2*pi
            inf = mod(intVal.inf, 2*pi);
            sup = mod(intVal.sup, 2*pi);

            %inf in [0, pi]
            if inf <= pi
                if sup < inf %due to mod computation  %I changed <= to < to solve the sin[0, 0] = (-1, 1 ) problem. 06.01.2016 (DG)
                    res.inf = -1;
                    res.sup = 1;
                elseif sup <= pi
                    res.inf = cos(sup);
                    res.sup = cos(inf);
                else 
                    res.inf = -1;
                    res.sup = max(cos(inf),cos(sup));
                end 
            %inf in [pi, 2*pi]
            else
                if sup <= pi
                    res.inf = min(cos(inf),cos(sup));
                    res.sup = 1;
                elseif sup < inf %due to mod computation
                    res.inf = -1;
                    res.sup = 1;
                else 
                    res.inf = cos(inf);
                    res.sup = cos(sup);
                end 
            end
        end

else

% matrix case
% rand(100, 100)
%time_CORA =
%   0.003018115402849
%time_INTLAB =
%   0.016251994873230

    % to preserve the shape    
    res = intVal;
    
    % find indices
    ind1 = find((intVal.sup - intVal.inf) >= 2*pi);   
    res.inf(ind1) = -1;
    res.sup(ind1) = 1;
    
    %remove multiples of 2*pi
    inf = mod(intVal.inf, 2*pi);
    sup = mod(intVal.sup, 2*pi);
    
    % inf in [0, pi]
    
    ind2 = find(((intVal.sup - intVal.inf) < 2*pi) & inf <= pi & sup < inf);    
    res.inf(ind2) = -1;
    res.sup(ind2) = 1;
    
    ind3 = find(((intVal.sup - intVal.inf) < 2*pi) & inf <= pi & sup <= pi & sup >= inf);    
    res.inf(ind3) = cos(sup(ind3));
    res.sup(ind3) = cos(inf(ind3));
    
    ind4 = find(((intVal.sup - intVal.inf) < 2*pi) & inf <= pi & sup > pi);
    res.inf(ind4) = -1;
    res.sup(ind4) = max(cos(inf(ind4)), cos(sup(ind4)));
    
    % inf in [pi, 2pi]
    
    ind5 = find(((intVal.sup - intVal.inf) < 2*pi) & inf > pi & sup > pi & sup < inf);    
    res.inf(ind5) = -1;
    res.sup(ind5) = 1;
    
    ind6 = find(((intVal.sup - intVal.inf) < 2*pi) & inf > pi & sup <= pi );
    res.inf(ind6) = min(cos(inf(ind6)), cos(sup(ind6)));
    res.sup(ind6) = 1;
    
    ind7 = find(((intVal.sup - intVal.inf) < 2*pi) & inf > pi & sup > pi & sup >= inf);    
    res.inf(ind7) = cos(inf(ind7));
    res.sup(ind7) = cos(sup(ind7));

       
end

%------------- END OF CODE --------------