function res = zonoBundle(obj)
% conZonotope - convert a mptPolytope object into a zonotope bundle object
%
% Syntax:  
%    res = zonoBundle(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    res - zonoBundle object
%
% Example: 
%    A = [-1 0;0 -1;1 1];
%    b = [1;1;1];
%    poly = mptPolytope(A,b);
%    zB = zonoBundle(poly);
%
%    figure
%    hold on
%    plot(poly,[1,2],'b','Filled',true,'EdgeColor','none');
%    xlim([-1.5,2.5]);
%    ylim([-1.5,2.5]);
%
%    figure
%    hold on
%    plot(zB,[1,2],'r','Filled',true,'EdgeColor','none');
%    xlim([-1.5,2.5]);
%    ylim([-1.5,2.5]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonoBundle/zonoBundle

% Author:       Niklas Kochdumper
% Written:      4-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % get halfspaces and vertices
    V = obj.P.V';
    A = obj.P.A;
    
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
    res = zonoBundle(Z);

%------------- END OF CODE --------------