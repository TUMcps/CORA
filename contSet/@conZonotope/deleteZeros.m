function res = deleteZeros(obj)
% deleteZeros - delete all-zero constraints and generators
%
% Syntax:  
%    res = deleteZeros(obj)
%
% Inputs:
%    obj - conZonotope object
%
% Outputs:
%    res - resulting conZonotope object
%
% Example: 
%    Z = [0 1 0 0 1;0 1 0 2 -1];
%    A = [-2 0 1 -1; 0 0 0 0];
%    b = [2;0];
%    cZ = conZonotope(Z,A,b);
%    
%    cZ_ = deleteZeros(cZ);
%
%    figure; hold on;
%    plot(cZ,[1,2],'r');
%    
%    figure; hold on;
%    plot(cZ_,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/deleteZeros
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author:       Niklas Kochdumper
% Written:      04-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % get object properties
    c = obj.Z(:,1); G = obj.Z(:,2:end); A = obj.A; b = obj.b;

    % remove all-zero constraints
    if ~isempty(A)
       ind = find(sum(abs([A,b]),2) < eps);
       A(ind,:) = []; b(ind) = [];
    end

    % remove all-zero generators
    ind = find(sum(abs([G;A]),1) < eps);
    G(:,ind) = []; A(:,ind) = [];
    
    % construct resulting conZonotope object
    res = conZonotope(c,G,A,b);        
    
%------------- END OF CODE --------------