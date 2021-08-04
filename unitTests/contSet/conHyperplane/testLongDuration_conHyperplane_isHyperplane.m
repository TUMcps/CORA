function res = testLongDuration_conHyperplane_isHyperplane
% testLongDuration_conHyperplane_isHyperplane - unit test function of isHyperplane
%
% Syntax:  
%    res = testLongDuration_conHyperplane_isHyperplane
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      22-March-2021
% Last update:  --- 
% Last revision:---

%------------- BEGIN CODE --------------
res = true;
nTests = 100;
for i=1:nTests
    n = randi(30);
    % empty constraint set
    hyp1 = conHyperplane(randn(1,n),randn,zeros(1,n),abs(randn));
    % constraint set not empty but still unconstrained hyperplane
    I = eye(n);
    a = abs(randn);
    hyp2 = conHyperplane(I(1,:),0,[I(1,:);-I(1,:)],[a;a]);
    if ~isHyperplane(hyp1) || ~isHyperplane(hyp2)
        res = false;
        break;
    end
    % constraint set which actually constrains hyperplane
    % could be 1D and then not consistent (see "conHyperplane"
    % constructor); otherwise trivially always an unconstrained hyperplane
    if n>1
        d = randn(1,n); d = d/norm(d);
        hyp3 = conHyperplane(d,1,[I;-I],ones(2*n,1));
        if isHyperplane(hyp3)
            res = false;
            break;
        end
    end
end

if res
    disp('testLongDuration_conHyperplane_isHyperplane successful');
else
    disp('testLongDuration_conHyperplane_isHyperplane failed');
end

%------------- END OF CODE --------------