function v_new = vertices_(obj,varargin)
% vertices_ - computes the vertices of a mptPolytope
%
% Syntax:  
%    V = vertices_(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    V - vertices object
%
% Example: 
%    C = [1 0 -1 0 1; 0 1 0 -1 1]';
%    d = [3; 2; 3; 2; 1];   
%    poly = mptPolytope(C,d);
%
%    V = vertices(poly)
%
%    figure; hold on;
%    plot(poly);
%    plot(V(1,:),V(2,:),'.k','MarkerSize',20);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/vertices

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      01-February-2011
% Last update:  12-June-2015
%               12-August-2016
%               09-May-2018 (NK, changed algorithm for vertex computation)
%               29-June-2018 (MA, all polytopes considered)
% Last revision:27-March-2023 (MW, rename vertices_)

%------------- BEGIN CODE --------------

% check if vertex representation already exists
if obj.P.hasVRep()
    
    v_new = obj.P.V;
    v_new = v_new';
    
else

    obj.P.minHRep();

    % for each polytope
    v = []; %init set of vertices
    for i=1:length(obj.P)
        try
            v_new = lcon2vert(obj.P(i).A,obj.P(i).b,[],[],[],0);
        catch ME
            % check for subspaces
            [~,S] = isFullDim(obj);
    
            % is polytope just a point?
            if isempty(S)
                v_new = center(obj)';
            elseif size(S,2) < size(obj.P.A,2)
                % polytope is non-degenerate in subspace;
                % we plug in
                %    x = S*y + c in Ax <= b, resulting in ASy <= b - Ac,
                % where y is of the subspace dimension and c is any point
                % within the original polytope P (we use the center-function)
                c = center(obj);
                P_ = mptPolytope(obj.P.A*S,obj.P.b - obj.P.A*c);
                % compute vertices for y in subspace
                v_new = vertices_(P_,'lcon2vert');
                % map back via x = S*y + c
                v_new = S*v_new + c;
                v_new = v_new';
                return
            else
                rethrow(ME);
            end
        end
        v = [v; v_new];
    end
    v_new = v';
end

%------------- END OF CODE --------------