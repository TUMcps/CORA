function S = and(hyp,S)
% and - computes the intersection of a constrained hyperplane with a set
%
% Syntax:  
%    S = and(hyp,S)
%
% Inputs:
%    hyp - conHyperplane object
%    S - contSet object
%
% Outputs:
%    S - contSet object
%
% Example: 
%    P = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    hyp = conHyperplane([1 1],2,[-1 0],-1);
%
%    res = hyp & P;
%
%    figure; hold on; xlim([-2,4]); ylim([-4,4]);
%    plot(hyp,[1,2],'r','LineWidth',3);
%    plot(P,[1,2],'b');
%    plot(res,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/and

% Author:       Niklas Kochdumper
% Written:      26-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% convert second set to polytope, otherwise infinite loop
if isa(S,'conHyperplane')
    S = mptPolytope(S);
end

% input argument check happens there
S = S & hyp;

%------------- END OF CODE --------------