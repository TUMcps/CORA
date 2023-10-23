function res = optBernstein( obj )
% optBernstein - range bound using bernstein polynomials
%
% Syntax:
%    res = optBernstein( obj )
%
% Inputs:
%    obj - taylm
%
% Outputs:
%    res - interval
%
% Example:
%   f = @(x) 1 + x.^5 - x.^4;
%   x = interval(0,1);
%   t = taylm(x,10,'x');
%   T = f(t);
%
%   intReal = interval(0.91808,1)
%
%   int = interval(T)
%   intBnb = interval(T,'bnb')
%   intBnbAdv = interval(T,'bnbAdv'), 
%   intLinQuad = interval(T,'linQuad')
%   intBernstein = interval(T, 'bernstein')
%
%   X = 0:0.01:1;
%   Y = f(X);
%   plot(X,Y);
%   xlim([0,1]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: taylm, taylm/interval, optBnbAdv, optLinQuad

% Authors:       Niklas Kochdumper
% Written:       03-February-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
   
    % convert polynomial part of the taylor model to a bernstein polynomial
    p = size(obj.monomials,2)-1;
    dom = interval(-ones(p,1),ones(p,1));

    B = poly2bernstein(obj.coefficients',obj.monomials(:,2:end)',dom);

    % compute enclosing interval
    res = interval(min(min(B)),max(max(B))) + obj.remainder;

% ------------------------------ END OF CODE ------------------------------
