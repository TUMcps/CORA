function [res,cert,scaling] = priv_venumZonotope(E,Z,tol,scalingToggle)
% priv_venumZonotope - checks whether an ellipsoid contains a
%    zonotope
%
% Syntax:
%    E = priv_venumZonotope(E,Z)
%
% Inputs:
%    E - ellipsoid object (circumbody)
%    Z - zonotope object (inbody)
%    tol - tolerance
%    scalingToggle - boolean allowing to choose whether scaling should be
%                    computed (which takes more time in general)
%
% Outputs:
%    res - true/false
%    cert - returns true iff the result of res could be
%           verified. For example, if res=false and cert=true, Z is
%           guaranteed to not be contained in E, whereas if res=false and
%           cert=false, nothing can be deduced (Z could still be
%           contained in E).
%           If res=true, then cert=true.
%    scaling - returns the smallest number 'scaling', such that
%           scaling*(E - E.center) + E.center contains Z.
%
% References:
%    [1] Kulmburg A, Brkan I, Althoff A, "Search-based and Stochastic
%        Solutions to the Zonotope and Ellipsotope Containment Problems",
%        ECC 2024
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/contains_

% Authors:       Adrian Kulmburg
% Written:       21-January-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

scaling = 0;
cert = true; % Since venumSearch is exact, cert is always true

% Let c1, c2 be the centers of Z, E. We prepare the norm-function,
% returning the norm of v-c2 w.r.t. the E-norm, where v is a given vertex
% of Z. Since v = c1 +- g_1 +- ... +- g_m, where g_i are the generators of
% Z, the norm of v-c2 is the same as the norm of G*nu + c1-c2, where
% G is the generator matrix of Z, nu = [+-1;...;+-1]. However,
% here, we strategically do NOT include the center just yet.
norm_E = @(p) ellipsoidNorm(E, p);

G = Z.generators;
% Instead, we add the combined centers to the generators of the inbody:
G = [Z.center-E.center G];

G_size = size(G,2);

% We now compute the norm of each vector in G, and sort the vectors in
% descending order; We keep the computed values for later use
generator_norms = zeros([1 G_size]);
for i = 1:G_size
   generator_norms(i) = norm_E(G(:,i));
end
[generator_norms,indices] = sort(generator_norms, 'descend');
G = G(:, indices);

% We also setup a heuristic based on the triangle inequality, that allows
% us to reject certain nodes that cannot possibly lead to a good result
heuristic = @(nu, value) value + generator_norms* (abs(abs(nu) - 1));
% Let me explain this line: it looks at all the coordinates of nu that are
% zero, and for those it adds the corresponding values of generator_norms.
% Then, it adds to that 'value', which is the current value of the node. By
% the triangle inequality, all nodes in the subtree of the node in question
% must have a value that is <= heuristic, which may be used to remove
% certain nodes of whose subtree cannot possibly lead to a high enough
% value to disprove containment.

% We now setup the queue of nodes to go through. We basically need to check
% all points of the form nu_1*h_1 + ... + nu_m*h_m, where h_i are the
% columns of H and nu_i are +1 or -1. These choices generate a binary tree.
% We then have a queue to save the nu's we still have to check
queue_nu = [];
% We also need to save the values of the nodes, which we put into a
% separate queue (but the two queues should be seen as linked)
queue_values = [];

% We start at the top:
queue_nu = [queue_nu zeros([G_size 1])];
queue_values = [queue_values 0];

% We can now go ahead with the iteration:
%   While we still have to check a point, keep the iteration going:
while not(isempty(queue_nu))
    % pop the last values of the queue:
    current_nu = queue_nu(:,end);         queue_nu(:,end) = [];
    current_value = queue_values(end);   queue_values(end) = [];
    
    % find first zero of current nu
    i = find(current_nu == 0, 1, 'first');
    if isempty(i)
        % If every component of current_nu is nonzero, we have hit a leaf
        % of the tree, and so we can stop there
        continue
    end
    % We now compute the children of the current_nu node
    child_positive = current_nu; child_negative = current_nu;
    child_positive(i) = 1;
    child_negative(i) = -1;
    
    % We compute the values of both children...
    child_positive_value = norm_E(G*child_positive);
    child_negative_value = norm_E(G*child_negative);
    
    % ... and order them, before adding them to the queue
    if child_positive_value <= child_negative_value
        queue_nu = [queue_nu child_positive child_negative];
        queue_values = [queue_values child_positive_value child_negative_value];
    else
        queue_nu = [queue_nu child_negative child_positive];
        queue_values = [queue_values child_negative_value child_positive_value];
    end
    
    % If any of the children has a value that is >1+tol, we have disproven
    % containment. By the ordering above, we just need to check the last
    % element of the queue
    if ~scalingToggle && queue_values(end) > 1+tol
        res = false;
        return
    else
        scaling = max([scaling queue_values(end)]);
    end
    
    % We now need to perform some 'cleaning' of the queue. First of all, we
    % need to compute the heuristic value of the sub-trees of both children
    % This part changes a bit depending on whether we need the scaling
    % factor, or whether we are only interested in a containment check
    if scalingToggle
        % If any of the two (or both) have a heuristic <= scaling (which is
        % the currently maximal value we have found so far), we don't need
        % to investigate the subtrees, they have no chance of beating the
        % maximum.
        if heuristic(queue_nu(:,end-1), queue_values(end-1)) <= scaling
            % So if the 'weakest' of the two nodes can be deleted from the
            % queue, do that
            queue_nu(:,end-1) = [];
            queue_values(end-1) = [];
            % Now, the question is whether the 'strongest' of the two also
            % needs to be deleted
            if heuristic(queue_nu(:,end), queue_values(end)) <= scaling
                queue_nu(:,end) = [];
                queue_values(end) = [];
            end

            % If something has been deleted, we don't need to sort anything in
            % the next step
            continue;
        end
        % No other situation can happen; it cannot happen, that the 'strongest'
        % node can be deleted, while the 'weakest' one doesn't; this follows
        % from the definition of the heuristic
    else
        % If any of the two (or both) have a heuristic <= 1+tol, there is no
        % need to keep them
        if heuristic(queue_nu(:,end-1), queue_values(end-1)) <= 1+tol
            % So if the 'weakest' of the two nodes can be deleted from the
            % queue, do that
            queue_nu(:,end-1) = [];
            queue_values(end-1) = [];
            % Now, the question is whether the 'strongest' of the two also
            % needs to be deleted
            if heuristic(queue_nu(:,end), queue_values(end)) <= 1+tol
                queue_nu(:,end) = [];
                queue_values(end) = [];
            end

            % If something has been deleted, we don't need to sort anything in
            % the next step
            continue;
        end
        % No other situation can happen; it cannot happen, that the 'strongest'
        % node can be deleted, while the 'weakest' one doesn't; this follows
        % from the definition of the heuristic
    end
    
    % Finally, we need to sort the 'weakest' of the two children, so that
    % the queue remains sorted (by the triangle inequality, the 'strongest'
    % child always has a node value that is at least as large as that of
    % the parent, so it may remain the last element of the queue)
    j = find(queue_values(1:end-1) >= queue_values(end-1), 1, 'first');
    queue_nu = [queue_nu(:,1:j-1) queue_nu(:,end-1) queue_nu(:,j:end-2) queue_nu(:,end)];
    
end

% If we could not find a node that disproves containment, there isn't any;
% so, the zonotope is contained in the other one
if scalingToggle
    res = scaling <= 1+tol;
else
    res = true;
end
end

% ------------------------------ END OF CODE ------------------------------
