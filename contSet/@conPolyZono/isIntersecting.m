function res = isIntersecting(cPZ,S,varargin)
% isIntersecting - determines if a constrained polynomial zonotope
%    intersects a set
%
% Syntax:  
%    res = isIntersecting(cPZ,S)
%    res = isIntersecting(cPZ,S,type)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    c = [0;0];
%    G = [8 4; 0 8];
%    expMat = [1 0; 0 1; 0 0];
%    A = [1 1 -0.25];
%    b = 0.75;
%    expMat_ = [2 0 0; 0 2 0; 0 0 1];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_);
%    pZ1 = polyZonotope([0;0],0.5*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%    pZ2 = polyZonotope([0;0],[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
%
%    res1 = isIntersecting(cPZ,pZ1,'approx')
%    res2 = isIntersecting(cPZ,pZ2,'approx')
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','b','Splits',15);
%    plot(pZ1,[1,2],'g');
%    plot(pZ2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/isIntersecting, isempty

% Author:       Niklas Kochdumper
% Written:      04-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[resFound,vars] = pre_isIntersecting('conPolyZono',cPZ,S,varargin{:});

% check premature exit
if resFound
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    % assign values
    cPZ = vars{1}; S = vars{2}; type = vars{3};
end


if strcmp(type,'exact')
    throw(CORAerror('CORA:noExactAlg',cPZ,S));
end

% % % get polyZonotope object
% % if ~isa(cPZ,'conPolyZono')
% %     temp = cPZ; cPZ = S; S = temp;
% % end

% call function for other set representations
if isa(S,'conPolyZono') || isa(S,'polyZonotope') || ...
   isa(S,'capsule')
    
    % fast check based on zonotope enclosure
    res = isIntersecting(zonotope(cPZ),zonotope(S));
    
    if ~res
        return; 
    end
    
    % convert second set to constrained polynomial zonotope
    S = conPolyZono(S);
    
    % compute intersection of the two sets
    int = cPZ & S;
    
    % check if the intersection is empty
    res = ~isempty(int);
    
    
elseif isa(S,'halfspace') || isa(S,'conHyperplane') || ...
       isa(S,'mptPolytope') || isa(S,'zonotope') || ...
       isa(S,'interval') || isa(S,'zonoBundle') || ...
       isa(S,'conZonotope') || isa(S,'ellipsoid')

    res = isIntersecting(S,cPZ,type);

else
    
    % throw error for given arguments
    throw(CORAerror('CORA:noops',cPZ,S));
end

%------------- END OF CODE --------------