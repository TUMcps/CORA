function res = cos(I)
% cos - Overloaded 'cos()' operator for intervals
%
% inf is x infimum, sup is x supremum
%
% [-1, 1]                       if (sup - inf) >= 2*pi,
% [-1, 1]                       if (sup - inf) < 2*pi and inf <= pi and sup <= pi and sup < inf,
% [cos(sup), cos(inf)]          if (sup - inf) < 2*pi and inf <= pi and sup <= pi and sup >= inf,
% [-1, max(cos(inf),cos(sup))]  if (sup - inf) < 2*pi and inf <= pi and sup > pi,
% [-1, 1]                       if (sup - inf) < 2*pi and inf > pi and sup > pi and sup < inf,
% [min(cos(inf),cos(sup)), 1]   if (sup - inf) < 2*pi and inf > pi and sup <= pi,
% [cos(inf), cos(sup)]          if (sup - inf) < 2*pi and inf > pi and sup > pi and sup >= inf.
%
% Syntax:
%    res = cos(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - interval object
%
% Example:
%    I = interval([-2;3],[3;4]);
%    res = cos(I);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff, Dmitry Grebenyuk
% Written:       25-June-2015
% Last update:   06-January-2016 (DG)
%                05-February-2016 (MA)
%                22-February-2016 (DG, the matrix case is rewritten)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% scalar case
if isnumeric(I)
    
    res = interval();

    %sup - inf >= 2pi
    if I.sup - I.inf >= 2*pi
        res.inf = -1;
        res.sup = 1;
    else
        %remove multiples of 2*pi
        lb = mod(I.inf, 2*pi);
        ub = mod(I.sup, 2*pi);

        %inf in [0, pi]
        if lb <= pi
            if ub < lb
                %due to mod computation
                %I changed <= to < to solve the sin[0, 0] = (-1, 1 ) problem. 06.01.2016 (DG)
                res.inf = -1;
                res.sup = 1;
            elseif ub <= pi
                res.inf = cos(ub);
                res.sup = cos(lb);
            else 
                res.inf = -1;
                res.sup = max(cos(lb),cos(ub));
            end 
        %inf in [pi, 2*pi]
        else
            if ub <= pi
                res.inf = min(cos(lb),cos(ub));
                res.sup = 1;
            elseif ub < lb %due to mod computation
                res.inf = -1;
                res.sup = 1;
            else 
                res.inf = cos(lb);
                res.sup = cos(ub);
            end 
        end
    end

else

    % to preserve the shape    
    res = I;
    
    % find indices
    ind1 = find((I.sup - I.inf) >= 2*pi);   
    res.inf(ind1) = -1;
    res.sup(ind1) = 1;
    
    %remove multiples of 2*pi
    lb = mod(I.inf, 2*pi);
    ub = mod(I.sup, 2*pi);
    
    % inf in [0, pi]
    
    ind2 = find(((I.sup - I.inf) < 2*pi) & lb <= pi & ub < lb);    
    res.inf(ind2) = -1;
    res.sup(ind2) = 1;
    
    ind3 = find(((I.sup - I.inf) < 2*pi) & lb <= pi & ub <= pi & ub >= lb);    
    res.inf(ind3) = cos(ub(ind3));
    res.sup(ind3) = cos(lb(ind3));
    
    ind4 = find(((I.sup - I.inf) < 2*pi) & lb <= pi & ub > pi);
    res.inf(ind4) = -1;
    res.sup(ind4) = max(cos(lb(ind4)), cos(ub(ind4)));
    
    % inf in [pi, 2pi]
    
    ind5 = find(((I.sup - I.inf) < 2*pi) & lb > pi & ub > pi & ub < lb);    
    res.inf(ind5) = -1;
    res.sup(ind5) = 1;
    
    ind6 = find(((I.sup - I.inf) < 2*pi) & lb > pi & ub <= pi );
    res.inf(ind6) = min(cos(lb(ind6)), cos(ub(ind6)));
    res.sup(ind6) = 1;
    
    ind7 = find(((I.sup - I.inf) < 2*pi) & lb > pi & ub > pi & ub >= lb);    
    res.inf(ind7) = cos(lb(ind7));
    res.sup(ind7) = cos(ub(ind7));
       
end

% ------------------------------ END OF CODE ------------------------------
