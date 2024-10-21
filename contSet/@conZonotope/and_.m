function res = and_(cZ,S,varargin)
% and_ - computes the intersection of a constrained zonotope with
%    other set representations
%
% Syntax:
%    res = and_(cZ,S)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object
%
% Outputs:
%    res - conZonotope object
%
% Example: 
%    % constrained zonotopes
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ2 = conZonotope(Z,A,b);
%
%    % halfspace and constrained hyperplane
%    P1 = polytope([1,-2],1);
%    P2 = polytope([-2 -0.5;1 0],[-4.25;2.5],[1,-2],1);
%
%    % compute intersection
%    res1 = cZ1 & cZ2;
%    res2 = cZ2 & P1;
%    res3 = cZ1 & P2;
%
%    % visualization
%    figure; hold on;
%    plot(cZ1,[1,2],'r');
%    plot(cZ2,[1,2],'b');
%    plot(res1,[1,2],'FaceColor','g');
%    title('Constrained zonotope');
%
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(P1,[1,2],'r','FaceAlpha',0.5);
%    plot(res2,[1,2],'FaceColor','g');
%    plot(cZ2,[1,2],'b');
%    title('halfspace');
%
%    figure; hold on; xlim([0,4]); ylim([-3,4]);
%    plot(P2,[1,2],'g');
%    plot(cZ1,[1,2],'r');
%    plot(res3,[1,2],'b','LineWidth',2);
%    title('Constrained hyperplane'); 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:       13-May-2018
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: 27-March-2023 (MW, rename and_)
%                28-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% call function with lower precedence
if isa(S,'contSet') && S.precedence < cZ.precedence
    res = and_(S,cZ,varargin{:});
    return
end

% constrained zonotope
if isa(S,'conZonotope')
    res = aux_and_conZonotope(cZ,S);
    return
end

% higher precedence cases
if isa(S,'zonotope') || isa(S,'interval') || isa(S,'zonoBundle')
    % convert to constrained zonotope
    res = aux_and_conZonotope(cZ,conZonotope(S));
    return
end

end


% Auxiliary functions -----------------------------------------------------

function res = aux_and_conZonotope(cZ,S)

% Calculate intersection according to equation (13) at Proposition 1 in
% reference paper [1]
Z = [cZ.c, cZ.G, zeros(size(S.G))];
A = blkdiag(cZ.A,S.A);
A = [A; cZ.G, -S.G];
b = [cZ.b; S.b; S.c - cZ.c];

res = conZonotope(Z,A,b);

% delete all zero constraints and generators
res = compact_(res,'zeros',eps);

end

% ------------------------------ END OF CODE ------------------------------
