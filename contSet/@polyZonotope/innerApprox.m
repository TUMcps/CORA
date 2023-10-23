function S = innerApprox(pZ,varargin)
% innerApprox - returns an inner-approximation of a polynomial zonotope
%    represented by a union of zonotopes
%
% Syntax:
%    S = innerApprox(pZ)
%    S = innerApprox(pZ,tol)
%
% Inputs:
%    pZ - polyZonotope object
%    tol - tolerance
%
% Outputs:
%    L - cell-array storing the zonotopes whose union inner-approximates
%        the polynomial zonotope
%
% Example:
%   pZ = polyZonotope([0;0],[2 0 2;0 2 2],[0.5;0],[1 0 3;0 1 1]);
%
%   S = innerApprox(pZ,1e-2);
%
%   figure; hold on;
%   plot(pZ,[1,2],'FaceColor','r');
%   for i = 1:length(S)
%       plot(S{i}); 
%   end
%
% References: 
%   [1] A. Goldsztejn and L. Jaulin. "Inner approximation of the range of 
%       vector-valued functions", 2010.
%   [2] O. Mullier and et al. "General Inner Approximation of Vector-valued 
%       Functions", Reliable Computing, 2013
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: innerApproxImage

% Authors:       Niklas Kochdumper
% Written:       21-December-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % compute inner-approximation of the dependent part using the 
    % algorithms in [1] and [2]
    p = size(pZ.E,1);
    f = @(x) aux_funcPoly(x,pZ);
    D = interval(-ones(p,1),ones(p,1));
    
    S = innerApproxImage(f,D,varargin{:});
    
    % add independent part
    for i = 1:length(S)
        if isempty(pZ.GI)
            S{i} = zonotope(S{i}); 
        else
            S{i} = zonotope(S{i}) + zonotope(zeros(dim(pZ),1),pZ.GI); 
        end
    end

end


% Auxiliary functions -----------------------------------------------------

function val = aux_funcPoly(x,pZ)
% nonlinear function that defines the dependent part of polynomial zonotope

    % initialization
    val = pZ.c;
    
    % dependent generators
    for i = 1:size(pZ.G,2)
        val = val + pZ.G(:,i) * prod(x.^pZ.E(:,i)); 
    end
end

% ------------------------------ END OF CODE ------------------------------
