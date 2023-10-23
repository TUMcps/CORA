function res = isIntersecting_(ls,S,type,varargin)
% isIntersecting_ - determines if a level set intersects a set using
%    the method described in Sec. 4.1 in [1]
%
% Syntax:
%    res = isIntersecting_(ls,S)
%    res = isIntersecting_(ls,S,type)
%
% Inputs:
%    ls - levelSet object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    syms x y
%    eq = sin(x) + y;
%    ls = levelSet(eq,[x;y],'<=');
%
%    I1 = interval([0.7;-0.3],[1.3;0.3]);
%    I2 = interval([0.2;-0.3],[0.8;0.3]);
%
%    isIntersecting(ls,I1,'approx')
%    isIntersecting(ls,I2,'approx')
%
%    figure; hold on; xlim([-1.5,1.5]); ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(I1,[1,2],'FaceColor','r');
%
%    figure; hold on; xlim([-1.5,1.5]); ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(I2,[1,2],'FaceColor','g');
%
% References:
%   [1] N. Kochdumper et al. "Reachability Analysis for Hybrid Systems with 
%       Nonlinear Guard Sets", HSCC 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, lfspace/isIntersecting_, conHyperplane/isIntersecting_

% Authors:       Niklas Kochdumper
% Written:       19-July-2019
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% currently, no exact method
if strcmp(type,'exact')
    throw(CORAerror('CORA:noops',ls,S));
end

% interval over-approximation
I = interval(S);

% evaluate non-linear function with interval arithmetic
I = ls.funHan(I);

% read infimum and supremum
ub = supremum(I);
lb = infimum(I);

% multiple inequality constraints or not
if ~iscell(ls.compOp)

    % switch case for the different types of level sets
    if strcmp(ls.compOp,'==')
        res = (ub > 0 | withinTol(ub,0)) & (lb < 0 | withinTol(lb,0)); 
    elseif strcmp(ls.compOp,'<=')
         res = lb < 0 | withinTol(lb,0);
    else
        res = lb < 0;
    end

else

    resVec = false(length(ls.compOp),1);

    % loop over all inequality constraints
    for i = 1:length(ls.compOp)

        if strcmp(ls.compOp{i},'<=')
            resVec(i) = lb(i) < 0 | withinTol(lb(i),0);
        else
            resVec(i) = lb(i) < 0;
        end
    end

    res = all(resVec);
end

% ------------------------------ END OF CODE ------------------------------
