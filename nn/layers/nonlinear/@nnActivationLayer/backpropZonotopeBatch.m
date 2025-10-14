function [gc, gG] = backpropZonotopeBatch(obj, c, G, gc, gG, options, updateWeights)
% backpropZonotopeBatch - compute the backpropagation for the previous input
%    with batches of zonotopes
%
% Syntax:
%    [gc,gG] = layeri.backpropZonotopeBatch(c,G,gc,gG,options,updateWeights);
%
% Inputs:
%    c, G - batch of input zonotopes; [n,q+1,b] = size([gc gG]),
%    gc, gG - batch of zonotope gradients; [n,q+1,b] = size([gc gG]),
%       where n is the number of dims, q the number of generators, and b the batch size
%    options - training parameters
%    updateWeights - only relevent for layer with learnable parameters
%
% Outputs:
%    gc, gG - zonotope gradients w.r.t the input
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork/backpropZonotopeBatch

% Authors:       Lukas Koller
% Written:       12-August-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Obtain indices of active generators.
genIds = obj.genIds;
% Obtain stored slope.
m = obj.backprop.store.coeffs;
% Obtain the number of dimensions and batch size.
[nk,bSz] = size(m);

% Check if there is an exact backpropagation through the image enclosure.
if options.nn.train.exact_backprop

    if ~options.nn.use_approx_error
        % TODO.
        throw(CORAerror('CORA:noExactAlg', obj, ...
            ['nnActivationLayer/backpropZonotopeBatch for ' ...
                'options.nn.train.exact_backprop && ' ...
                    '~options.nn.use_approx_error']))
    end

    % Obtain the stored gradients.
    m_l = obj.backprop.store.m_l; % grad of m w.r.t to l
    m_u = obj.backprop.store.m_u; % grad of m w.r.t to u
    % Obtain the stored gradients.
    dl_l = obj.backprop.store.dl_l; % grad of dl w.r.t to l
    dl_u = obj.backprop.store.dl_u; % grad of dl w.r.t to u
    du_l = obj.backprop.store.du_l; % grad of du w.r.t to l
    du_u = obj.backprop.store.du_u; % grad of du w.r.t to u

    % During the forward propagation we store the gradient of the slope m
    % and the approximation errors [dl,du] w.r.t. to the bounds l and u. 
    % From, u = c + r, l = c - r, and r = sum(abs(G),2) we can compute the 
    % gradients w.r.t. the centers and the generators.

    % The gradient of the radius is the sign of the generator entries, 
    % i.e., d(r)/d(G) = sign(G). Therefore, we compute the sign of the 
    % active generators.
    r_G = sign(G(:,genIds,:));

    % Now, we can compute the gradient of the slope w.r.t. the generators.
    % I.e., d(m)/d(G) = d(m)/d(l) * d(l)/d(G) + d(m)/d(u) * d(u)/d(G),
    % where d(l)/d(G) = -d(r)/d(G) and d(u)/d(G) = d(r)/d(G).
    m_G = permute(m_u - m_l,[1 3 2]).*r_G;

    % Obtain stored indices of the dimensions of the approximation 
    % errors.
    dDimsIdx = obj.backprop.store.dDimsIdx;
    notdDimsIdx = obj.backprop.store.notdDimsIdx;

    % Compute the ids of the dimensions for approximation errors.
    [dDims,~] = ind2sub([nk bSz],dDimsIdx);
    % Obtain the indices of the generators for the approximation errors.
    approxErrGenIds = obj.approxErrGenIds;
    % Obtain the number of considered approximation errors.
    dn = length(approxErrGenIds);

    if dn > 0
        % Compute indices for approximation errors into the gradients of 
        % the generator matrix.
        GdIdx = reshape(sub2ind(size(gG), ...
           reshape(dDims,1,[]),...
           repmat(approxErrGenIds,1,bSz), ...
           repelem(1:bSz,1,dn)),[dn bSz]);
    else
        % There are no approximation errors stored in the generator matrix.
        GdIdx = zeros([0 bSz],'like',gG);
    end

    if options.nn.store_approx_error_gradients
        % Compute the gradients w.r.t. the slope.
        % Obtain the gradient of the approximation error w.r.t the slope.
        dl_m = obj.backprop.store.dl_m;
        du_m = obj.backprop.store.du_m;
        % Compute the gradient of the output center w.r.t. the slope 
        % (derived from d(c')/d(m) = c + 1/2*(d(du)/d(m) + d(dl)/d(m)).
        % Compute the gradient of the offset.
        offset_m = 1/2*(du_m + dl_m);
        if options.nn.interval_center
            % Obtain the lower and upper bounds of the center.
            cl = reshape(c(:,1,:),[nk bSz]);
            cu = reshape(c(:,1,:),[nk bSz]);
            % The gradient of the offset is 0 for all approximation error 
            % that are stored in the interval center.
            offset_m(notdDimsIdx) = 0;
            % Set the gradient for the approximation errors stored in the
            % interval center.
            dcl_m = dl_m;
            dcl_m(dDimsIdx) = 0;
            dcu_m = du_m;
            dcu_m(dDimsIdx) = 0;
            % Compute the gradient.
            c_m = permute(cat(3, ...
                cl + offset_m - dcl_m, ...
                cu + offset_m + dcu_m),[1 3 2]);
        else
            % Compute the gradient.
            c_m = c + offset_m;
        end
        % Compute the gradient of the output generators w.r.t. the slope 
        % (derived from d(G')/d(m) = G, and for approximation error 
        % generators: d(G')/d(m) = 1/2*(d(du)/d(m) - d(dl)/d(m)).
        G_m = zeros(size(gG),'like',gG);
        G_m(:,genIds,:) = G(:,genIds,:);
        G_m(GdIdx) = 1/2*(du_m(dDimsIdx) - dl_m(dDimsIdx));
        % Compute the slope gradient from all generators.
        G_slopeGrad = reshape(sum(gG.*G_m,2),[nk bSz]);
        % Compute the slope gradient from the center.
        c_slopeGrad = gc.*c_m;
        if options.nn.interval_center
            % Sum the gradients of the bounds of the interval center.
            c_slopeGrad = reshape(sum(c_slopeGrad,2),[nk bSz]);
        end
        % Compute the slope gradient.
        slopeGrad = c_slopeGrad + G_slopeGrad;
        % Update the slope gradient.
        obj.backprop.store.slope_gradients = slopeGrad;

        % Add the slope gradient for optimization.
        if ~isfield(obj.backprop.store,'dm') ...
                || ~all(size(m) == size(obj.backprop.store.dm))
            obj.backprop.store.dm = zeros(size(slopeGrad),'like',slopeGrad);
        end
        % Obtain the step size.
        stepSize = options.nn.verif_slope_optim_step_size;
        % Compute the gradient step.
        % dmStep = stepSize*min(max(slopeGrad,0),1);
        dmStep = stepSize*(slopeGrad - min(slopeGrad,[],1))...
            ./(max(slopeGrad,[],1) - min(slopeGrad,[],1) + 1e-6);
        % Update the slope offset.
        obj.backprop.store.dm = obj.backprop.store.dm - dmStep;
    end

    % Permute the dimensions of the generator gradients for easier
    % indexing (switch generator dimension with input dimension).
    r_G = permute(r_G,[2 1 3]);

    % Compute gradients of the approximation errors w.r.t. the 
    % generators (analogous to the slope gradient computation).
    dl_G = reshape(dl_u(dDimsIdx) - dl_l(dDimsIdx),1,[]).*r_G(:,dDimsIdx);
    du_G = reshape(du_u(dDimsIdx) - du_l(dDimsIdx),1,[]).*r_G(:,dDimsIdx);

    % Compute outer product of the input and gradient generator matrix. 
    outProdG = sum(G(:,genIds,:).*gG(:,genIds,:),2);

    if options.nn.interval_center
        % Obtain the bounds of the gradient interval center.
        gl = gc(:,1,:);
        gu = gc(:,2,:);

        % For the interval center, l = cl - r and u = cu + r. Thus, 
        % d(l)/d(cl) = 1, d(u)/d(cl) = 0, d(l)/d(cu) = 0, and 
        % d(u)/d(cu) = 1.

        % Compute the gradient of the slope w.r.t. the interval center
        % bounds.
        m_cl = m_l;
        m_cu = m_u;

        % Compute gradients of the approximation errors w.r.t. the interval 
        % center (analogous to the slope gradient computation).
        dl_cl = dl_l(dDimsIdx);
        du_cl = du_l(dDimsIdx);
        dl_cu = dl_u(dDimsIdx);
        du_cu = du_u(dDimsIdx);

        % Compute the outer product of input and gradient centers.
        outProd_cl = c(:,1,:).*gl;
        outProd_cu = c(:,2,:).*gu;
        % Add the outer products.
        outProd = outProd_cl + outProd_cu + outProdG;

        % Compute the factors for the approximation errors.
        fdl_l = gl(dDimsIdx) - gG(GdIdx);
        fdu_l = gl(dDimsIdx) + gG(GdIdx);
        fdl_u = gu(dDimsIdx) - gG(GdIdx);
        fdu_u = gu(dDimsIdx) + gG(GdIdx);
        % Add factor for lower and upper approximation errors.
        fdl = fdl_u + fdl_l;
        fdu = fdu_u + fdu_l;

        % Compute the gradients w.r.t. the input center bounds (without 
        % the gradient for the approximation errors).
        gcl = m.*gl(:,:) + outProd(:,:).*m_cl;
        gcu = m.*gu(:,:) + outProd(:,:).*m_cu;
        % Add the gradients for the offset of through the approximation
        % errors that are added as generators.
        gcl(dDimsIdx) = gcl(dDimsIdx) + 1/2*fdu.*du_cl + 1/2*fdl.*dl_cl;
        gcu(dDimsIdx) = gcu(dDimsIdx) + 1/2*fdu.*du_cu + 1/2*fdl.*dl_cu;
        % Add gradient for the approximation errors that were added to the
        % interval center.
        gcl(notdDimsIdx) = gcl(notdDimsIdx) - dl_l(notdDimsIdx);
        gcu(notdDimsIdx) = gcu(notdDimsIdx) + du_u(notdDimsIdx);
        % Concatenate the center gradient bounds.
        gc = permute(cat(3,gcl,gcu),[1 3 2]);
    else
        % Compute the gradient of the slope w.r.t. the center (derived from
        % d(m)/d(c) = d(m)/d(l) * d(l)/d(c) + d(m)/d(u) * d(u)/d(c), where
        % d(l)/d(c) = d(u)/d(c) = 1;
        m_c = m_u + m_l;

        % Compute gradients of the approximation errors w.r.t. the center 
        % (analogous to the slope gradient computation).
        dl_c = dl_u(dDimsIdx) + dl_l(dDimsIdx);
        du_c = du_u(dDimsIdx) + du_l(dDimsIdx);

        % Compute the outer product of input and gradient center.
        outProdc = permute(c.*gc,[1 3 2]);
        % Add the outer products.
        outProd = outProdc + outProdG;

        % Compute the factors for the approximation errors.
        fdl = gc(dDimsIdx) - gG(GdIdx);
        fdu = gc(dDimsIdx) + gG(GdIdx);

        % Compute the gradients w.r.t. the input center.
        gc = m.*gc + outProd(:,:).*m_c;
        % Add the gradient for the approximation errors.
        gc(dDimsIdx) = gc(dDimsIdx) + 1/2*(fdu.*du_c + fdl.*dl_c);
    end

    if options.nn.store_approx_error_gradients
        % Extract the gradients.
        dgrad = zeros([nk bSz],'like',gG);
        dgrad(dDimsIdx) = gG(GdIdx);
        % Store the gradients w.r.t. the approximation errors.
        obj.backprop.store.approx_error_gradients = dgrad;
    end

    % Compute the gradients w.r.t. the input generators.
    gG(:,genIds,:) = permute(m,[1 3 2]).*gG(:,genIds,:) + outProd.*m_G;

    % Permute the dimensions of the generator gradients for easier
    % indexing (switch generator dimension with input dimension).
    gG_ = permute(gG(:,genIds,:),[2 1 3]);

    % Add the gradient for the approximation errors.
    gG_(:,dDimsIdx) = gG_(:,dDimsIdx) + 1/2*(fdu(:)'.*du_G + fdl(:)'.*dl_G);

    % Re-permute the dimensions of the generator gradients.
    gG(:,genIds,:) = permute(gG_,[2 1 3]);
else
    % Consider the linear approximation as fixed. Use the slope of the
    % approximation for backpropagation.
    if options.nn.interval_center
        gc = permute(m,[1 3 2]).*gc;
    else
        gc = gc.*m;
    end
    gG(:,genIds,:) = gG(:,genIds,:).*permute(m,[1 3 2]);
end
% Only keep the gradient with corresponding inputs.
gG(:,(max(genIds)+1):end,:) = [];

end

% ------------------------------ END OF CODE ------------------------------
