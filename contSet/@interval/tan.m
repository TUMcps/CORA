function I = tan(I)
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
%    I - interval object
%
% Example: 
%    I = interval(0.1,0.2);
%    tan(I)
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

% Authors:       Daniel Althoff, Dmitry Grebenyuk, Matthias Althoff, Mark Wetzlinger
% Written:       03-November-2015
% Last update:   14-January-2016 (DG)
%                05-February-2016 (MA)
%                17-March-2016 (DA, Speed improvement)
%                12-September-2016 (DG, fixed tan([-pi/2, ...]) )
% Last revision: 18-January-2024 (MW, faster algorithm)

% ------------------------------ BEGIN CODE -------------------------------

% copy to avoid constructor call
lb = I.inf;
ub = I.sup;

% init interval with -Inf/Inf
[n,m] = size(I);
I.inf = -Inf(n,m);
I.sup = Inf(n,m);

% only dimensions with a diameter larger smaller than pi or where
% tan(inf) < tan(sup) have non-Inf values
taninf = tan(lb);
tansup = tan(ub);
ind = (ub - lb) < pi & taninf <= tansup;
I.inf(ind) = taninf(ind);
I.sup(ind) = tansup(ind);


%% previous algorithm [Eq. (14), 1]
% res = I;   % to preserve the shape of the imput matrix
% 
% ind1 = I.sup - I.inf >= pi;   % find all intervals broader then pi
% res.inf(ind1) = -Inf;
% res.sup(ind1) = +Inf;
% 
% inf = mod(I.inf, pi);  % project the left intervals onto [0, pi]
% sup = mod(I.sup, pi);
% 
% ind1 = inf < pi/2 & (sup < inf | sup > pi/2) & (I.sup - I.inf < pi);
% res.inf(ind1) = -Inf;
% res.sup(ind1) = +Inf;
%                     
% ind1 = inf < pi/2 & (sup >= inf & sup <= pi/2) & (I.sup - I.inf < pi);
% res.inf(ind1) = tan(inf(ind1));
% res.sup(ind1) = tan(sup(ind1));
% 
% ind1 = inf >= pi/2 & (sup < inf & sup > pi/2) & (I.sup - I.inf < pi);
% res.inf(ind1) = -Inf;
% res.sup(ind1) = +Inf;
%                     
% ind1 = inf >= pi/2 & (sup >= inf | sup <= pi/2) & (I.sup - I.inf < pi);
% res.inf(ind1) = tan(inf(ind1));
% res.sup(ind1) = tan(sup(ind1));
% 
% % a fix for a case when inf is -pi/2. In that case, the mod function produce
% % pi/2, which has an ambiguity of +Inf or -Inf at that point
% ind1 = res.sup < res.inf;
% res.inf(ind1) = -res.inf(ind1);

% ------------------------------ END OF CODE ------------------------------
