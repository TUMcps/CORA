function zB = zonoBundle(P)
% zonoBundle - convert a polytope object into a zonotope bundle object
%
% Syntax:
%    res = zonoBundle(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    zB - zonoBundle object
%
% Example: 
%    A = [-1 0;0 -1;1 1];
%    b = [1;1;1];
%    P = polytope(A,b);
%    zB = zonoBundle(P);
%
%    figure; hold on
%    xlim([-1.5,2.5]); ylim([-1.5,2.5]);
%    plot(P,[1,2],'FaceColor',colorblind('b'));
%
%    figure; hold on
%    xlim([-1.5,2.5]); ylim([-1.5,2.5]);
%    plot(zB,[1,2],'FaceColor',colorblind('r'));
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonoBundle/zonoBundle

% Authors:       Niklas Kochdumper
% Written:       04-November-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
% get halfspaces and vertices
V = vertices(P);
A = P.A;

% loop over all halfspaces
Z = cell(size(A,1),1);

for i = 1:size(A,1)
   
    % get basis orthogonal to halfspace
    B = gramSchmidt(A(i,:)');
    
    % transform vertices into new space
    V_ = B' * V;
    
    % compute enclosing box
    box = zonotope(interval(min(V_,[],2),max(V_,[],2)));
    
    % transform back to original space
    Z{i} = B*box;
end

% construct resulting zonoBundle object
zB = zonoBundle(Z);

% ------------------------------ END OF CODE ------------------------------
