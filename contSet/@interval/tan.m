function res = tan(I)
% tan - Overloaded 'tan()' operator for intervals
%
% inf is infimum, sup is supremum
%
% [-Inf, Inf]                   if (sup - inf) >= pi,
% [-Inf, Inf]                   if (sup - inf) < pi) and inf < pi/2 and (sup < inf or sup > pi/2),
% [tan(inf), tan(sup)]          if (sup - inf) < pi and inf < pi/2 and (sup >= inf and sup <= pi/2),
% [-Inf, Inf]                   if (sup - inf) < pi) and inf >= pi/2 and (sup < inf and sup > pi/2),
% [tan(inf), tan(sup)]          if (sup - inf) < pi and inf >= pi/2 and (sup >= inf or sup <= pi/2).
%
% Syntax:  
%    res = tan(I)
%
% Inputs:
%    I - interval object
%
% Outputs:
%    res - interval object
%
% Example: 
%    I = interval(0.1,0.2);
%    tan(I)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Author:       Daniel Althoff, Dmitry Grebenyuk, Matthias Althoff
% Written:      03-November-2015
% Last update:  14-January-2016 (DG)
%               05-February-2016 (MA)
%               17-March-2016 (DA, Speed improvement)
%               12-September-2016 (DG, fixed tan([-pi/2, ...]) )
% Last revision:---

%------------- BEGIN CODE --------------

res = I;   % to preserve the shape of the imput matrix

ind1 = I.sup - I.inf >= pi;   % find all intervals broader then pi
res.inf(ind1) = -Inf;
res.sup(ind1) = +Inf;

inf = mod(I.inf, pi);  % project the left intervals onto [0, pi]
sup = mod(I.sup, pi);

ind1 = inf < pi/2 & (sup < inf | sup > pi/2) & (I.sup - I.inf < pi);
res.inf(ind1) = -Inf;
res.sup(ind1) = +Inf;
                    
ind1 = inf < pi/2 & (sup >= inf & sup <= pi/2) & (I.sup - I.inf < pi);
res.inf(ind1) = tan(inf(ind1));
res.sup(ind1) = tan(sup(ind1));

ind1 = inf >= pi/2 & (sup < inf & sup > pi/2) & (I.sup - I.inf < pi);
res.inf(ind1) = -Inf;
res.sup(ind1) = +Inf;
                    
ind1 = inf >= pi/2 & (sup >= inf | sup <= pi/2) & (I.sup - I.inf < pi);
res.inf(ind1) = tan(inf(ind1));
res.sup(ind1) = tan(sup(ind1));

% a fix for a case when inf is -pi/2. In that case, the mod function produce
% pi/2, which has an ambiguity of +Inf or -Inf at that point
ind1 = res.sup < res.inf;
res.inf(ind1) = -res.inf(ind1);
 
%------------- END OF CODE --------------