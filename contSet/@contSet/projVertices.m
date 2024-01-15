function V = projVertices(S,varargin)
% projVertices - computes the vertices of a 2D projection of a set based on
%    support function evaluation of the projected set; if no more support
%    vectors can be found, the algorithm terminates
%    this function also supports degenerate sets (lines, points)
%
% Syntax:
%    V = projVertices(S)
%
% Inputs:
%    S - contSet object
%    dims - dimensions for projection
%
% Outputs:
%    V - list of vertices in the projected space
%
% Example:
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ = conZonotope(Z,A,b);
%
%    V = projVertices(cZ);
%
%    figure; hold on;
%    plot(cZ);
%    scatter(V(1,:),V(2,:),16,'r','filled');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       21-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% set default values
dims = setDefaultValues({[1,2]},varargin);

% check input arguments (only supported for convex sets with edges)
inputArgsCheck({ {S,'att',{'conZonotope','interval','zonoBundle','zonotope'}} ...
                 {dims,'att','numeric',{'integer','positive','size',[1 2]}} });

% project the set
S = project(S,dims);

% init vertices
V = [];
if representsa_(S,'emptySet',1e-10)
    % no vertices if set is empty
    return
end

% other options for support function evaluation
otherOptions = {};
if isa(S,'polyZonotope') || isa(S,'conPolyZono')
    otherOptions = {'interval',8,1e-3};
end

% figure; blobsize = 24; subplot(1,2,1); hold on;

% compute support vectors of three directions

% 0 degrees
[~,V(:,1)] = supportFunc_(S,[1;0],'upper',otherOptions{:});
% plot([0,1],[0,0],'k');

% 120 degrees
angle_pi = 120*pi / 180;
dir = [cos(angle_pi) -sin(angle_pi); sin(angle_pi) cos(angle_pi)] * [1;0];
% plot([0,dir(1)],[0,dir(2)],'k');
[~,V(:,2)] = supportFunc_(S,dir,'upper',otherOptions{:});

% 240 degrees
angle_pi = 240*pi / 180;
dir = [cos(angle_pi) -sin(angle_pi); sin(angle_pi) cos(angle_pi)] * [1;0];
% plot([0,dir(1)],[0,dir(2)],'k');
[~,V(:,3)] = supportFunc_(S,dir,'upper',otherOptions{:});

% copy last point for easier indexing (will be deleted later)
V(:,4) = V(:,1);
% indices of second and third vertices
idx2 = 2; idx3 = 3;

% ensure that there are no duplicates (other than 1 = end)
if compareMatrices(V(:,3),V(:,1:2),1e-14,'subset')
    % vertex (240 degrees) is duplicate
    V(:,3) = []; idx3 = [];
end
if withinTol(V(:,2),V(:,1))
    V(:,2) = []; idx2 = []; idx3 = idx3-1;
end

% all sections between already computed vertices have to be investigated
sections = mat2cell([(1:size(V,2)-1)', (2:size(V,2))'],ones(size(V,2)-1,1),2);

% split angles between neighboring directions until no new information
while ~isempty(sections)
    
%     subplot(1,2,2); scatter(V(1,:),V(2,:),blobsize,'k','filled'); hold on;

    % analyze first section in the list of sections
    section = sections{1};
    sections = sections(2:end);

    % direction = normal vector on halfspace spanned by points on section
    v = V(:,section(2)) - V(:,section(1));
    dir = [v(2); -v(1)];
    dir = dir / vecnorm(dir);

%     subplot(1,2,1); hold on; plot([0,dir(1)],[0,dir(2)],'r');
    
    % compute support vector for the new direction
    [~,V_new] = supportFunc_(S,dir,'upper',otherOptions{:});

%     subplot(1,2,2); scatter(V_new(1),V_new(2),blobsize,'r'); hold off;
    
    % compute vectors:
    % - from start vertex to computed vertex
    % - from computed vertex to end vertex
    ptsStartMidEnd = [V_new-V(:,section(1)), V(:,section(2))-V_new];

    % check whether sections is completed
    if compareMatrices(V_new,V,1e-6,'subset') ...
            || rank(ptsStartMidEnd,1e-6) < 2
        % new vertex is on a line with start and end points of the current
        % section -> discard vertex, section completed
%         subplot(1,2,1); plot([0,dir(1)],[0,dir(2)],'w');
    else
        % vertices are not on a line
%         subplot(1,2,1); plot([0,dir(1)],[0,dir(2)],'k');
        
        % add vertex to list at index between points from current section
        V = [V(:,1:section(1)) V_new V(:,section(2):end)];

        % shift indices of other sections
        for s=1:length(sections)
            sections{s} = sections{s} + 1;
        end
        % add two new sections at the beginning
        sections = [{[section(1) section(1)+1];[section(2) section(2)+1]}; sections];

        % shift indices of second and third vertex
        if section(2) <= idx2
            idx2 = idx2 + 1;
        end
        if section(2) <= idx3
            idx3 = idx3 + 1;
        end
    end

end

% remove last vertex (was only added for convenience)
V = V(:,1:end-1);

% check if first three vertices are actual vertices
V = aux_checkActualVertex(V,1);
V = aux_checkActualVertex(V,idx2);
V = aux_checkActualVertex(V,idx3);

% close;

end


% Auxiliary functions -----------------------------------------------------

function V = aux_checkActualVertex(V,idx)
% check if vertex at position idx is an actual vertex or if its on the line
% segment between the neighboring vertices

% idx is empty if it has been eliminated at the start
if isempty(idx) || size(V,2) <= 2
    return
end

% compute vectors to neighboing point in list
if idx == 1
    ptsStartMidEnd = [V(:,idx)-V(:,idx+1), V(:,end)-V(:,idx)];
elseif idx == size(V,2)
    ptsStartMidEnd = [V(:,idx)-V(:,1), V(:,idx-1)-V(:,idx)];
else
    ptsStartMidEnd = [V(:,idx)-V(:,idx+1), V(:,idx-1)-V(:,idx)];
end

% check rank
if rank(ptsStartMidEnd,1e-12) < 2
    V(:,idx) = [];
end

end

% ------------------------------ END OF CODE ------------------------------
