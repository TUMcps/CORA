function res = contractForwardBackward(f,dom)
% contractForwardBackward - implementation of the forward-backward
%                           contractor acc. to Sec. 4.2.4 in [1]
%
% Syntax:
%    res = contractForwardBackward(f,dom)
%
% Inputs:
%    f - function handle for the constraint f(x) = 0
%    dom - initial domain (class: interval)
%
% Outputs:
%    res - contracted domain (class: interval)
%
% Example: 
%    f = @(x) x(1)^2 + x(2)^2 - 4;
%    dom = interval([1;1],[3;3]);
%   
%    res = contract(f,dom,'forwardBackward');
%
%    figure; hold on
%    xlim([-3,3]); ylim([-3,3]);
%    plot(dom,[1,2],'r');
%    plot(res,[1,2],'g');
%    syms x1 x2
%    ls = levelSet(f([x1;x2]),[x1;x2],'==');
%    plot(ls,[1,2],'b');   
%
% References:
%    [1] L. Jaulin et al. "Applied Interval Analysis", 2006
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contract

% Authors:       Zhuoling Li, Niklas Kochdumper
% Written:       04-November-2019 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % compute syntax tree nodes for all variables
    vars = [];
    
    for i = 1:length(dom)
        vars = [vars;syntaxTree(dom(i),i)];
    end

    % forward iteration: compute syntax tree
    synTree = f(vars);
    
    % backward iteration
    res = dom;
    
    for i = 1:size(synTree,1)
        try 
            res = backpropagation(synTree(i),interval(0,0),res);
        catch ex
            if strcmp(ex.identifier,'CORA:emptySet')
                res = [];
                return
            else
                rethrow(ex);
            end
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
