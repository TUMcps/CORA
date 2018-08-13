function [bound,xMin,domMin,xMax,domMax] = globVerBounds(func,dom,tol,varargin)
% globVerBounds - determine the gloabl bounds (minimum and maximum) of a
%                 function on the given search domain up to a certain
%                 precision
%
% Syntax:  
%    [val,xOpt,domOpt] = globVerBounds(func,dom,tol)
%    [val,xOpt,domOpt] = globVerBounds(func,dom,tol,max_order,opt_method,eps,tolerance)
%
% Inputs:
%    func - function for which the bounds are computed (provided as a 
%           function handle)
%    dom - multi-dimensional search domain (class: interval)
%    tol - tolerance for the determined values of the bounds. The real
%          minimum "min_real" is located inside the following interval:
%          min_real \in [val, val + tol]. The same logic applies for the
%          maximum
%    max_order - the maximal polynomial order of a monomial stored in a 
%                polynomial part of a taylor model
%    opt_method - method used to calculate interval over-approximations of
%                 taylor models during the calculation of the initial
%                 taylor model object
%                  'int': standard interval arithmetic (default)
%                  'bnb': branch and bound method is used to find min/max
%                  'bnbAdv': branch and bound with re-expansion of taylor models
%                  'linQuad': optimization with Linear Dominated Bounder (LDB)
%                             and Quadratic Fast Bounder (QFB)
%   eps - precision for the selected optimization method for the talyor 
%         over-approximation with intervals (opt_method = 'bnb', 
%         opt_method = 'bnbAdv' and opt_method = 'linQuad')
%   tolerance - taylor model monomials with coefficients smaller than this
%               value are moved to the remainder
%
% Outputs:
%    bound - determined bounds of the function on the search domain
%            (class: interval)
%    xMin - best guess for the point at which the function reaches it's
%           minimum
%    domMin - domain in which the determined lower bound of the function 
%             minimum is located  
%    xMax - best guess for the point at which the function reaches it's
%           maximum
%    domMax - domain in which the determined upper bound of the function 
%             maximum is located  
%
% Example:
%   % 1D-function with maximum 1 and minimum 0.91808 at x=0.8
%   f= @(x) 1 + x.^5 - x.^4;
%   x = interval(0,1);
%   [int,xMin,domMin,xMax,domMax] = globVerBounds(f,x,0.0001);
%
%   % plot the functions and the calculated boundaries
%   X = infimum(x):0.01:supremum(x);
%   Y = f(X);
%   plot(X,Y,'b');
%   hold on;
%   plot(X,ones(size(X))*infimum(int),'r');
%   plot(X,ones(size(X))*supremum(int),'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: globVerMinimization, optLinQuad
%
% References: 
%   [1] K. Makino et al. "Verified Global Optimization with Taylor Model 
%       based Range Bounders"
%   [2] K. Makino "Rigorous analysis of nonlinear motion in particle
%       accelerators"

% Author:       Niklas Kochdumper
% Written:      16-April-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % Implementation of the Verified Global Optimizer concept from
    % reference paper [1]

    % default values
    max_order = 10;
    opt_method = 'int';
    eps = 0.001;
    tolerance = 1e-8;
    
    % parse input arguments
    if nargin < 3
       error('Wrong syntax!. Type "help globVerMinimization" for help.'); 
    end
    if nargin >= 4 && ~isempty(varargin{1})
       max_order = varargin{1}; 
    end
    if nargin >= 5 && ~isempty(varargin{2})
       opt_method = varargin{2}; 
    end
    if nargin >= 6 && ~isempty(varargin{3})
       eps = varargin{3}; 
    end
    if nargin >= 7 && ~isempty(varargin{4})
       tolerance = varargin{4}; 
    end
    
    % create taylor models
    t = taylm(dom,max_order,'x',opt_method,eps,tolerance);
    T = func(t);
    
    % create affine arithmetic objects
    a = affine(dom,'x',opt_method,eps,tolerance);
    A = func(a);
    
    % Minimum
    [minVal,xMin,domMin] = globVerMinimization(func,dom,tol,max_order,opt_method,eps,tolerance,T,A);
    
    % Maximum
    func_ = @(x) -func(x);
    [maxVal,xMax,domMax] = globVerMinimization(func_,dom,tol,max_order,opt_method,eps,tolerance,-T,-A);
    
    % overall bounds
    bound = interval(minVal,-maxVal);
    
    %------------- END OF CODE --------------