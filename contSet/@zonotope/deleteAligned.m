function Z = deleteAligned(Z)
% deleteAligned - combines aligned generators to a single generator;
%    a tolerance is used to determine alignment, so this function
%    does not necessarily return an over-approximation of the original
%    zonotope--for this, use reduce instead
%
% Syntax:  
%    Z = deleteAligned(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    Z - zonotope object
%
% Example:
%    Z1 = zonotope([1;0],[1 0 1 1 -2 0; 0 1 0 1 -2 -3]);
%    Z2 = deleteAligned(Z1);
% 
%    plot(Z1); hold on;
%    plot(Z2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:        Matthias Althoff
% Written:       15-January-2009
% Last update:   27-Aug-2019
% Last revision: ---

%------------- BEGIN CODE --------------

%extract center and generator matrix
c=center(Z);
G=generators(Z);

%Delete zero-generators
G=nonzeroFilter(G);

%normalize generators
G_norm = G./vecnorm(G);

% tolerance for alignment
tol = 1-1e-3;

%find equal generators
i = 1;
while i < length(G(1,:))
    G_act = G_norm(:,i);
    ind = find(abs(G_act'*G_norm(:,(i+1):end)) > tol);
    if ~isempty(ind)
        ind = ind+i;
        for iAdd = 1:length(ind)
            %add generators
            G(:,i) = G(:,i) + sign(G_act'*G_norm(:,ind(iAdd)))*G(:,ind(iAdd));
        end
        for iAdd = 1:length(ind)
            %remove generators
            G(:,ind(iAdd)) = [];
            G_norm(:,ind(iAdd)) = [];
            %change ind to correct for cancellation
            ind = ind - 1;
        end
    end  
    %increase i
    i = i + 1;
end

Z.Z = [c,G];

%------------- END OF CODE --------------
