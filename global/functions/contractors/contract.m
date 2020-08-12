function res = contract(f,dom,varargin)
% contract - contracts a interval domain to tightly enclose nonlinear
%            constraints
%
% Syntax:  
%    res = contract(f,dom)
%    res = contract(f,dom,alg)
%
% Inputs:
%    f - function handle for the constraint f(x) = 0
%    dom - initial domain (class: interval)
%    alg - algorithm used for contraction. The available algorithms are 
%          'forwardBackward' (see Sec. 4.2.4 in [1]),  'linearize' 
%          (see Sec. 4.3.4 in [1]), and 'polyBox' (see [2])
%
% Outputs:
%    res - contracted domain (class: interval)
%
% Example: 
%    f = @(x) x(1)^2 + x(2)^2 - 4;
%    dom = interval([1;1],[3;3]);
%   
%    res = contract(f,dom);
%
%    figure
%    hold on
%    xlim([-3,3]);
%    ylim([-3,3]);
%    plot(dom,[1,2],'r');
%    plot(res,[1,2],'g');
%    syms x1 x2
%    ls = levelSet(f([x1;x2]),[x1;x2],'==');
%    plot(ls,[1,2],'b');   
%
% References:
%    [1] L. Jaulin et al. "Applied Interval Analysis", 2006
%    [2] G. Trombettoni et al. "A Box-Consistency Contractor Based on 
%        Extremal Functions", 2010
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contractForwardBackward, contractParallelLinearization

% Author:       Niklas Kochdumper
% Written:      04-November-2019 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    alg = 'forwardBackward';
    
    if nargin >= 3
       alg = varargin{1}; 
    end

    % contract the domain using the selected algorithm
    if strcmp(alg,'forwardBackward')
        res = contractForwardBackward(f,dom);
    elseif strcmp(alg,'linearize')
        res = contractParallelLinearization(f,dom);
    elseif strcmp(alg,'polyBox')
        res = contractPolyBoxRevise(f,dom);
    end
end

%------------- END OF CODE --------------