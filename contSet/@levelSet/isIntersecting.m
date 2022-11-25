function res = isIntersecting(obj1,obj2,varargin)
% isIntersecting - determines if levelSet obj1 intersects obj2 with the
%                  method described in Sec. 4.1 in [1]
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%    res = isIntersecting(obj1,obj2,type)
%
% Inputs:
%    obj1 - halfspace object
%    obj2 - conSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - 1/0 if set is intersecting, or not
%
% Example: 
%    syms x y
%    eq = sin(x) + y;
%    ls = levelSet(eq,[x;y],'<=');
%
%    int1 = interval([0.7;-0.3],[1.3;0.3]);
%    int2 = interval([0.2;-0.3],[0.8;0.3]);
%
%    isIntersecting(ls,int1,'approx')
%    isIntersecting(ls,int2,'approx')
%
%    figure
%    hold on
%    xlim([-1.5,1.5]);
%    ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(int1,[1,2],'r','Filled',true,'EdgeColor','none');
%
%    figure
%    hold on
%    xlim([-1.5,1.5]);
%    ylim([-1,1]);
%    plot(ls,[1,2],'b');
%    plot(int2,[1,2],'g','Filled',true,'EdgeColor','none');
%
% References:
%   [1] N. Kochdumper et al. "Reachability Analysis for Hybrid Systems with 
%       Nonlinear Guard Sets", HSCC 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace/isIntersecting, conHyperplane/isIntersecting

% Author:       Niklas Kochdumper
% Written:      19-July-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    type = 'exact';
    
    if nargin > 2 && ~isempty(varargin{1})
       type = varargin{1};
    end
    
    if strcmp(type,'exact')
       error('No exact algorithm implemented for this set representation!');
    end

    % interval over-approximation
    int = interval(obj2);

    % evaluate non-linear function with interval arithmetic
    intRes = obj1.funHan(int);

    sup = supremum(intRes);
    infi = infimum(intRes);

    % multiple inequality constraints or not
    if ~iscell(obj1.compOp)

        % switch case for the different types of level sets
        if strcmp(obj1.compOp,'==')
            res = (sup >= 0) & (infi <= 0); 
        elseif strcmp(obj1.compOp,'<=')
            res = (infi <= 0);
        else
            res = (infi < 0);
        end

    else

        resVec = zeros(length(obj1.compOp),1);

        % loop over all inequality constraints
        for i = 1:length(obj1.compOp)

            if strcmp(obj1.compOp{i},'<=')
                resVec(i) = (infi(i) <= 0);
            else
                resVec(i) = (infi(i) < 0);
            end
        end

        res = all(resVec == 1);
    end

%------------- END OF CODE --------------