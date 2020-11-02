function p = randPointExtreme(P)
% randPointExtreme - generates a random extreme point of a polytope
%
% Syntax:  
%    p = randPointExtreme(P)
%
% Inputs:
%    P - mptPolytope object
%
% Outputs:
%    p - random point in R^n
%
% Example: 
%    poly = mptPolytope.generateRandom(2);
%    p = randPointExtreme(poly);
%
%    figure; hold on;
%    plot(poly,[1,2],'r');
%    plot(p(1),p(2),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: randPoint

% Author:       Niklas Kochdumper
% Written:      30-October-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % check if vertex representation is available
    if P.P.hasVRep
        
       % select random vertex
       V = P.P.V;
       ind = randi([1,size(V,1)]);
       p = V(ind,:)';
       
    else
       
       % center polytope at origin
       c = center(P);
       P = P - c;
       
       % select random direction
       n = length(c);
       d = rand(n,1) - 0.5*ones(n,1);
       d = d./norm(d);
        
       % compute farest point in this direction that is still in polytope
       [~,x] = supportFunc(P,d);
       p = x + c;
       
    end
end

%------------- END OF CODE --------------