function res = sin(I)
% sin - Overloaded 'sin()' operator for intervals
%
% inf is x infimum, sup is x supremum
%
% [-1, 1]                       if (sup - inf) >= 2*pi,
% [-1, 1]                       if (sup - inf) < 2*pi) and inf <= pi/2 and sup < inf),
% [sin(inf), sin(sup)]          if (sup - inf) < 2*pi and inf <= pi/2 and sup <= pi/2 and sup >= inf,
% [min(sin(inf), sin(sup)), 1]  if (sup - inf) < 2*pi and inf <= pi/2 and sup > pi/2 and sup <= 3/2*pi,
% [-1, 1]                       if (sup - inf) < 2*pi and inf <= pi/2 and sup > 3/2*pi),
% [-1, 1]                       if (sup - inf) < 2*pi and inf > pi/2 and inf <= 3/2*pi and sup > pi/2 and sup < inf,
% [-1, max(sin(inf), sin(sup))] if (sup - inf) < 2*pi and inf > pi/2 and inf <= 3/2*pi and sup <= pi/2,
% [sin(sup), sin(inf)]          if (sup - inf) < 2*pi and inf > pi/2 and inf <= 3/2*pi and sup > pi/2 and sup <= 3/2*pi and sup >= inf,
% [-1, 1]                       if (sup - inf) < 2*pi and inf > 3/2*pi and inf <= 2*pi and sup > 3/2*pi and sup < inf,
% [sin(inf), sin(sup)]          if (sup - inf) < 2*pi and inf > 3/2*pi and inf <= 2*pi and sup <= pi/2,
% [min(sin(inf), sin(sup)) , 1] if (sup - inf) < 2*pi and inf > 3/2*pi and inf <= 2*pi and sup > pi/2 and sup <= 3/2*pi,
% [sin(inf), sin(sup)]          if (sup - inf) < 2*pi and inf > 3/2*pi and inf <= 2*pi and sup > 3/2*pi and sup >= inf.
%
% Syntax:
%    res = sin(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - interval object
%
% Example: 
%    I = interval([-2;-1],[3;4]);
%    sin(I)
%
% References:
%    [1] M. Althoff, D. Grebenyuk, "Implementation of Interval Arithmetic
%        in CORA 2016", ARCH'16.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff, Dmitry Grebenyuk, Mark Wetzlinger
% Written:       24-June-2015
% Last update:   13-January-2016 (DG)
%                05-February-2016 (MA)
%                06-February-2016 (DG)
%                10-February-2016 (MA)
%                22-February-2016 (DG, the matrix case is rewritten)
%                10-January-2024 (MW, fix condition for scalar case)
% Last revision: 18-January-2024 (MW, call fast algorithm in interval/cos)

% ------------------------------ BEGIN CODE -------------------------------

res = cos(I - pi/2);

%% previous algorithm [Eq. (12), 1]
% % scalar case
% if isscalar(I)
%     
%     % init (overwritten below)
%     res = interval.Inf(1);
% 
%     %sup - inf >= 2pi
%     if I.sup - I.inf >= 2*pi
%         res.inf = -1;
%         res.sup = 1;
%     else
%         %remove multiples of 2*pi
%         inf = mod(I.inf, 2*pi);
%         sup = mod(I.sup, 2*pi);
% 
%         %inf in [0, pi/2]
%         if inf <= pi/2
%             %due to mod computation
%             %I changed <= to < to solve the sin[0, 0] = (-1, 1 ) problem.
%             if sup < inf
%                 res.inf = -1;
%                 res.sup = 1;
%             elseif sup <= pi/2
%                 res.inf = sin(inf);
%                 res.sup = sin(sup);
%             elseif sup <= 3/2*pi
%                 res.inf = min(sin(inf), sin(sup));
%                 res.sup = 1;
%             else 
%                 res.inf = -1;
%                 res.sup = 1;
%             end
%         %inf in [pi/2, 3/2*pi]
%         elseif inf <= 3/2*pi
%             if sup <= pi/2
%                 res.inf = -1;
%                 res.sup = max(sin(inf), sin(sup));
%             elseif sup < inf %due to mod computation
%                 res.inf = -1;
%                 res.sup = 1;
%             elseif sup <= 3/2*pi
%                 res.inf = sin(sup);
%                 res.sup = sin(inf);
%             else 
%                 res.inf = -1;
%                 res.sup = max(sin(inf), sin(sup));
%             end    
%         %inf in [3/2*pi, 2*pi]
%         else
%             if sup <= pi/2
%                 res.inf = sin(inf);
%                 res.sup = sin(sup);
%             elseif sup <= 3/2*pi %this was wrong
%                 res.inf = min(sin(inf), sin(sup));
%                 res.sup = 1;
%             elseif sup < inf %this was wrong
%                 res.inf = -1;
%                 res.sup = 1;
%             else 
%                 res.inf = sin(inf);
%                 res.sup = sin(sup);
%             end 
%         end
%     end
%     
% else
% 
%     % to preserve the shape    
%     res = I;
%     
%     % find indices
%     ind1 = find((I.sup - I.inf) >= 2*pi);   
%     res.inf(ind1) = -1;
%     res.sup(ind1) = 1;
%     
%     %remove multiples of 2*pi
%     inf = mod(I.inf, 2*pi);
%     sup = mod(I.sup, 2*pi);
%     
%     % inf in [0, pi/2]
%     
%     ind2 = find(((I.sup - I.inf) < 2*pi) & inf <= pi/2 & sup < inf);    
%     res.inf(ind2) = -1;
%     res.sup(ind2) = 1;
%     
%     ind3 = find(((I.sup - I.inf) < 2*pi) & inf <= pi/2 & sup <= pi/2 & sup >= inf);    
%     res.inf(ind3) = sin(inf(ind3));
%     res.sup(ind3) = sin(sup(ind3));
%     
%     ind4 = find(((I.sup - I.inf) < 2*pi) & inf <= pi/2 & sup > pi/2 & sup <= 3/2*pi);
%     res.inf(ind4) = min(sin(inf(ind4)), sin(sup(ind4)));
%     res.sup(ind4) = 1;
%     
%     ind5 = find(((I.sup - I.inf) < 2*pi) & inf <= pi/2 & sup > 3/2*pi);
%     res.inf(ind5) = -1;
%     res.sup(ind5) = 1;
%     
%     % inf in [pi/2, 3/2*pi]
%     
%     ind6 = find(((I.sup - I.inf) < 2*pi) & inf > pi/2 & inf <= 3/2*pi & sup > pi/2 & sup < inf);    
%     res.inf(ind6) = -1;
%     res.sup(ind6) = 1;
%     
%     ind7 = find(((I.sup - I.inf) < 2*pi) & inf > pi/2 & inf <= 3/2*pi & sup <= pi/2 );
%     res.inf(ind7) = -1;
%     res.sup(ind7) = max(sin(inf(ind7)), sin(sup(ind7)));
%     
%     ind8 = find(((I.sup - I.inf) < 2*pi) & inf > pi/2 & inf <= 3/2*pi & sup > pi/2 & sup <= 3/2*pi & sup >= inf);    
%     res.inf(ind8) = sin(sup(ind8));
%     res.sup(ind8) = sin(inf(ind8));
%     
%     ind9 = find(((I.sup - I.inf) < 2*pi) & inf > pi/2 & inf <= 3/2*pi & sup > 3/2*pi & sup >= inf);    
%     res.inf(ind9) = -1;
%     res.sup(ind9) = max(sin(inf(ind9)), sin(sup(ind9)));
%     
%     % inf in [3/2*pi, 2*pi]
%     
%     ind10 = find(((I.sup - I.inf) < 2*pi) & inf > 3/2*pi & inf <= 2*pi & sup > 3/2*pi & sup < inf);    
%     res.inf(ind10) = -1;
%     res.sup(ind10) = 1;
%     
%     ind11 = find(((I.sup - I.inf) < 2*pi) & inf > 3/2*pi & inf <= 2*pi & sup <= pi/2);    
%     res.inf(ind11) = sin(inf(ind11));
%     res.sup(ind11) = sin(sup(ind11));
%     
%     ind12 = find(((I.sup - I.inf) < 2*pi) & inf > 3/2*pi & inf <= 2*pi & sup > pi/2 & sup <= 3/2*pi );
%     res.inf(ind12) = min(sin(inf(ind12)), sin(sup(ind12)));
%     res.sup(ind12) = 1;
%     
%     ind13 = find(((I.sup - I.inf) < 2*pi) & inf > 3/2*pi & inf <= 2*pi & sup > 3/2*pi & sup >= inf);    
%     res.inf(ind13) = sin(inf(ind13));
%     res.sup(ind13) = sin(sup(ind13));
%        
% end

% ------------------------------ END OF CODE ------------------------------
