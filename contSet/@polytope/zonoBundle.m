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
% Last update:   03-January-2023 (MW, handle unbounded case)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if representsa_(P,'fullspace',0)
    % conversion of fullspace object not possible
    throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
        'can therefore not be converted into a zonotope bundle.']));
end

% get halfspaces and vertices
try
    V = vertices(P);
catch ME
    if ~isempty(P.bounded.val) && ~P.bounded.val
        throw(CORAerror('CORA:specialError',['Polytope is unbounded and '...
            'can therefore not be converted into a zonotope bundle.']));
    end
    rethrow(ME);
end

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
