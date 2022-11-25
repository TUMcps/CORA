function res = test_conHyperplane_isHyperplane
% test_conHyperplane_isHyperplane - unit test function of isHyperplane
%
% Syntax:  
%    res = test_conHyperplane_isHyperplane
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

% Author:       Victor Gassmann, Mark Wetzlinger
% Written:      27-July-2021
% Last update:  --- 
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% empty constraint set
a = [3; 2; -1];
b = 0.5;
C = [0 0 0];
d = 1;
hyp1 = conHyperplane(a,b,C,d);

if ~isHyperplane(hyp1)
    res = false;
end

% constraint set not empty but still unconstrained hyperplane
I = eye(2);
a = 0.7;
hyp2 = conHyperplane(I(1,:),0,[I(1,:);-I(1,:)],[a;a]);

if ~isHyperplane(hyp2)
    res = false;
end
    
% constraint set which actually constrains hyperplane
% could be 1D and then not consistent (see "conHyperplane"
% constructor); otherwise trivially always an unconstrained hyperplane
d = [1; -1]; d = d/norm(d);
hyp3 = conHyperplane(d,1,[I(1,:);-I(1,:)],ones(2,1));

if isHyperplane(hyp3)
    res = false;
end

    
if res
    disp('test_conHyperplane_isHyperplane successful');
else
    disp('test_conHyperplane_isHyperplane failed');
end

%------------- END OF CODE --------------