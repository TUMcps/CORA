classdef nnGCNLayer < nnGNNLayer
% nnGCNLayer - class for graph convolutional layer (gcn)
%
% Syntax:
%    obj = nnGCNLayer(name)
%
% Inputs:
%    name - name of the layers
%
% Outputs:
%    obj - generated GCNlayer
%
% References:
%    [1] Kipf, T. N., et al. (2016). Semi-supervised classification with
%        graph convolutional networks. arXiv preprint arXiv:1609.02907.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Gerild Pjetri, Tianze Huang, Filip Dimovski, Tobias Ladner
% Written:       13-December-2022
% Last update:   23-February-2023
%                13-February-2024 (TL, moved W,b into nnGNNLinearLayer)
%                29-February-2024 (TL, saved MP, correct ids)
%                20-March-2024 (TL, faster matrix quad map implementation)
%                16-April-2024 (TL, polynomial matrix zonotope)
% Last revision: 15-January-2023
%                27-August-2023 (TL, practical course revisions)

% ------------------------------ BEGIN CODE -------------------------------

properties (Constant)
    is_refinable = false
    message_passing = setproperty()
end

methods
    % constructor
    function obj = nnGCNLayer(name)
        if nargin < 1
            name = [];
        end
        
        % call super class constructor
        obj@nnGNNLayer(name)
    end

    function outputSize = getOutputSize(obj, inputSize, G)
        outputSize = inputSize;
    end

    function [nin, nout] = getNumNeurons(obj)
        nin = []; nout = [];
    end
end

% evaluate ----------------------------------------------------------------

methods (Access = {?nnLayer, ?neuralNetwork})


    % numeric
    function r = evaluateNumeric(obj, input, options)
        % compute message passing
        MP = aux_computeMessagePassing(obj, size(input,1), options,1);
        
        % message passing
        r = vecleftmtimes(MP,input);
    end

    % sensitivity
    function S = evaluateSensitivity(obj, S, x, options)
        % compute message passing
        MP = aux_computeMessagePassing(obj, size(x,1), options);

        % init
        nrNodes = options.nn.graph.numnodes;
        nrFeat = size(x,1)/nrNodes;

        % message passing
        S = full(S * kron(speye(nrFeat),MP));
        
    end

    % zonotope/polyZonotope
    function [c, G, GI, E, id, id_, ind, ind_] = evaluatePolyZonotope(obj, c, G, GI, E, id, id_, ind, ind_, options)
        % message passing
        MP = aux_computeMessagePassing(obj, numel(c), options, id_);

        % init
        nrNodes = options.nn.graph.numnodes;
        nrFeat = numel(c)/nrNodes;

        % check if we have uncertainty in the message passing
        if isnumeric(MP)
            % graph with no perturbed edges
            c = vecleftmtimes(MP,c);
            G = vecleftmtimes(MP,G);
            GI = vecleftmtimes(MP,GI);

        else
            % graph with perturbed edges

            % create a polyZonotope from the sent arguments
            X = polyZonotope(c, G, GI, E, id);

            % calculate the output from mQuadMap(MP,X)
            Y = obj.aux_mQuadMap_gcn(MP, X, nrNodes, nrNodes, nrFeat);

            % extract properties
            c = Y.c;
            G = Y.G;
            GI = Y.GI;
            if isempty(GI)
                GI = zeros(length(c), 0);
            end
            E = Y.E;
            id = Y.id;
            id_ = max(id);
            ind = find(prod(ones(size(E))-mod(E, 2), 1) == 1);
            ind_ = setdiff(1:size(E, 2), ind);
        end

        % make it full again
        c = full(c);
        G = full(G);
        GI = full(GI);
    end

end


% Auxiliary functions -----------------------------------------------------

methods (Access = protected)

    function MP = aux_computeMessagePassing(obj, dim, options, maxId)
        % compute vectorized message passing
        G = options.nn.graph;

        % read uncertain edges
        idx_pert_edges = options.nn.idx_pert_edges;

        % check if MP is saved
        DAD_saved = nnGCNLayer.message_passing.val;

        % recompute MP
        if isempty(idx_pert_edges)
            % no edge is perturbed

            if isnumeric(DAD_saved) && ~isempty(DAD_saved)
                % reuse message passing
                DAD = DAD_saved;
            else
                % recompute DAD

                % obtain adjacency matrix
                A = G.adjacency;
                
                % compute \tilde{D}^(-1/2)
                % can be computed element-wise as D is diagonal
                D = diag(sqrt(1./sum(A)));
                if any(isinf(D),"all")
                    throw(CORAerror('CORA:outOfDomain','validDomain','D\\in(0,Inf]^{N \\times N}'));
                end
                D(isinf(D)) = 0;
    
                % compute D^(-1/2) * A * D^(-1/2)
                DAD = D * A * D;
            end

            % assign message passing
            MP = DAD;

        else
            % uncertain graph structure

            if isa(DAD_saved,'polyZonotope')
                % reuse message passing
                DAD = DAD_saved;
            else
                % recompute DAD

                % 1. get uncertain adjacency matrix
                numNodes = G.numnodes;
    
                % find perturbed edges
                [s, t] = findedge(G, idx_pert_edges);
                numPertEdges = numel(s);
    
                % perturb adjacency matrix
                A_G = zeros(numNodes^2*numPertEdges,1);
                A_G(numNodes^2*((1:numPertEdges)'-1) + numNodes*(t-1) + s) = 0.5;
                A_G(numNodes^2*((1:numPertEdges)'-1) + numNodes*(s-1) + t) = 0.5;
                A_G = reshape(A_G,numNodes^2,numPertEdges);
    
                % center is 0.5 for all perturbed edges 
                % and original adjacency matrix otherwise.
                A_c = sum(A_G,2);
                A_c = A_c + (A_c == 0) .* reshape(G.adjacency,[],1);
    
                % create polyZonotope of uncertain adjacency matrix
                A = polyZonotope(A_c,A_G);
    
                % 2. compute diagonal entries of degree matrix
                D_diag = kron(ones(1,numNodes),eye(numNodes)) * A;
    
                % 3. enclose inverse square root
                invSqrtLayer = nnInvSqrtRootLayer();
                invSqrtLayer.order = options.nn.invsqrt_order;
                options.nn.use_approx_error = options.nn.invsqrt_use_approx_error;
                D_diag_invsqrt = invSqrtLayer.evaluate(D_diag,options);
    
                % 4. compute full degree matrix
                In = speye(numNodes^2);
                idx = 1:(numNodes+1):numNodes^2;
                D_invsqrt = In(:,idx) * D_diag_invsqrt;
    
                % 5. extend exponent matrix of A
                p_diff = size(D_invsqrt.E,1) - size(A.E,1);
                A = polyZonotope(A.c,A.G,A.GI,[A.E;zeros(p_diff,size(A.E,1))]);
    
                % 6. compute message passing
                DA = obj.aux_mQuadMap_gcn(D_invsqrt,A,numNodes,numNodes,numNodes);
                DAD = obj.aux_mQuadMap_gcn(DA,D_invsqrt,numNodes,numNodes,numNodes);
    
                % fix ids
                DAD = DAD.replaceId(maxId + DAD.id);
            end

            % not required to adapt to match \vec{X}
            MP = DAD;

            % keyboard
        end

        DAD_saved = nnGCNLayer.message_passing;
        DAD_saved.val = DAD;
    end

     function PZ = aux_mQuadMap_gcn(obj, PZ1, PZ2, n, k, m)
        % compute quadratic map of uncertain matrices
        % correctly considers ids of given polynomial zonotopes

        % transform given sets to polynomial matrix zonotopes

        % init
        c1 = PZ1.c; G1 = PZ1.G; h1 = size(G1,2);
        c2 = PZ2.c; G2 = PZ2.G; h2 = size(G2,2);

        % reshape to respective matric
        c1 = reshape(c1,n,k); G1 = reshape(G1,n,k,h1,1);
        c2 = reshape(c2,k,m); G2 = reshape(G2,k,m,1,h2);

        % compute respective matrix multiplications
        c1c2 = c1 * c2;
        G1c2 = pagemtimes(G1, c2);
        c1G2 = pagemtimes(c1, G2);
        G1G2 = pagemtimes(G1, G2);

        % reshape back to vector
        c1c2 = reshape(c1c2,n*m,[]);
        G1c2 = reshape(G1c2,n*m,[]);
        c1G2 = reshape(c1G2,n*m,[]);
        G1G2 = reshape(G1G2,n*m,[]);

        % compute output generator matrix
        G = [G1c2, c1G2, G1G2];

        % compute output exponent matrix
        id1 = PZ1.id; id2 = PZ2.id;
        id = unique([id1;id2]);
        % extend exponent matrices
        E1 = (id1' == id) * PZ1.E;
        E2 = (id2' == id) * PZ2.E;
        E = [E1, E2];
        if ~isempty(PZ1.E) && ~isempty(PZ2.E)
            % sum each column of E1 to E2 and concatenate
            E = [E, reshape(E1 + reshape(E2,[],1,h2),[],h1*h2)];
        end

        % init resulting polynomial zonotope
        PZ = polyZonotope(c1c2, G, [], E, id);
    end

end

end

% ------------------------------ END OF CODE ------------------------------
