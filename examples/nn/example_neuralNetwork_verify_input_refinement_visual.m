function res = example_neuralNetwork_verify_input_refinement_visual()
% example_neuralNetwork_verify_input_refinement_visual - example for 
%   visualizing the iterative input refinement used in 
%   neuralNetwork/verify.
%
% Syntax:
%    res = example_neuralNetwork_verify_input_refinement_visual()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       20-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Reset the random number generator.
rng('default');

% Generate random neural network.
nn = neuralNetwork.generateRandom(NrInputs=2,NrOutputs=2, ...
    ActivationFun='relu',NrLayers=3,NrHiddenNeurons=4);
nn.layers(end) = [];

% Specify initial set.    
cx = [0; 0];
Gx = [1 0; 0 1];
Ax = zeros([0 2]);
bx = zeros([0 1]);
X = zonotope(cx,Gx);
% Obtain number of generators.
[~,q] = size(Gx);

% Compute the exact output set.
cYs = aux_splitSet(nn,cx,Gx,Ax,bx);
% Unify all sets.
Y = polygon(emptySet(2));
for i=1:length(cYs)
    Yi = aux_2ConZonoWithEqConst(cYs{i}.c,cYs{i}.G,cYs{i}.A,cYs{i}.b);
    if ~representsa(Yi,'emptySet')
        Y = Y | polygon(Yi);
    end
end

% Specify unsafe set.
% C = [-5 1; 1 1];       
% d = [-1.8; -1.2];  
C = [-1 1];    
d = [-2.3];
safeSet = false;
U = polytope(C,d);

% Specify number of refinement iterations.
nIter = 5;

% Initialize input set.
Xi = X;

% Initialize figure;
f = aux_initFigure(Y);
% Plot initial input set.
f = aux_plotSet(f,Xi,1);
% Plot specification.
f = aux_plotSet(f,U,2,CORAcolor('CORA:unsafe'));

for i=1:nIter
    % Enclose the output set.
    [cy,Gy,Ay,by] = aux_evaluateConZonotope(nn,cx,Gx,Ax,bx);

    % Construct the output set.
    Yi = aux_2ConZonoWithEqConst(cy,Gy,Ay,by);

    % Plot output set.
    f = aux_plotSet(f,Yi,2);
    
    % Compute intersection withe the unsafe set.
    Ay = [Ay; C*Gy];
    by = [by; d - C*cy];
   
    % Construct the intersection.
    uYi = aux_2ConZonoWithEqConst(cy,Gy,Ay,by);
    if representsa(uYi,'emptySet')
        f = aux_plotSet(f,Yi,2,CORAcolor('CORA:safe'));
        break;
    end

    % Plot initial input set.
    f = aux_plotSet(f,uYi,2,CORAcolor('CORA:purple'));

    % Add zero-generators.
    Gx = [Gx zeros([size(Gx,1) size(Gy,2) - size(Gx,2)])];
    Ax = [Ax zeros([size(Ax,1) size(Ay,2) - size(Ax,2)])];

    % Apply constraints to the to input set.
    Ax = [Ax; Ay];
    bx = [bx; by];

    % Construct new constraint zonotope.
    Xi = aux_2ConZonoWithEqConst(cx,Gx,Ax,bx);
    % Plot new input set.
    f = aux_plotSet(f,Xi,1);
end

end


% Auxiliary functions -----------------------------------------------------

function [c,G,A,b] = aux_evaluateConZonotope(nn,c,G,A,b)
    % Loop through the layers.
    for i=1:length(nn.layers)
        % Obtain the i-the layer.
        layeri = nn.layers{i};

        if isa(layeri,'nnLinearLayer')
            % Obtain the weight matrix and bias vector.
            Wi = layeri.W;
            bi = layeri.b;
            % Apply the affine map.
            c = Wi*c + bi;
            G = Wi*G;
            % Constraints are not changed.
        elseif isa(layeri,'nnReLULayer')
            % Obtain number of dimensions.
            [n,~] = size(G);
            % Obtain number of constraints.
            [p,~] = size(A);

            % Compute bounds of the input set.
            [l,u] = aux_boundsOfConZonotope(c,G,A,b,'exact');

            % Compute approximation slope.
            mi = (layeri.f(u) - layeri.f(l))./(u - l);
            % Compute approximation errors.
            xs = [l zeros(n,1) u];
            ys = layeri.f(xs);
            % Compute approximation error at candidates.
            ds = ys - mi.*xs;
            % Only consider extreme points within the bounds.
            notInBoundsIdx = (xs < l | xs > u);
            % Obtain lower approximation error.
            ds(notInBoundsIdx) = inf;
            dl = min(ds,[],2,'linear');
            % Obtain upper approximation error.
            ds(notInBoundsIdx) = -inf;
            du = max(ds,[],2,'linear');

            % Apply enclosure.
            c_ = mi.*c + 1/2*(du + dl);
            G_ = [mi.*G diag(1/2*(du - dl))];

            % Pad the constraints.
            A = [A zeros(p,n)];
            % Add constraints: 
            % (i) x >= 0 <--> -G'*beta - dr <= c'
            A = [A; -G_];
            b = [b; c_];
            % (ii) ReLU(x) >= x <--> (G - G')*beta - dr <= c' - c
            A = [A; [G zeros(n,n)] - G_];
            b = [b; c_ - c];

            % Update the center and generators.
            c = c_;
            G = G_;
        else
            throw(CORAerror('CORA:nnLayerNotSupported', ...
                layeri,'evaluateConZonotope'));
        end
    end
end

function [l,u] = aux_boundsOfConZonotope(c,G,A,b,varargin)
    % Set default bounding method.
    boundingMethod = setDefaultValues({'dual-iter'},varargin);

    % Construct the batched constrained zonotope.
    cZs = struct('c',c,'G',G,'dr',0,'A',A,'b',b);

    % Specify the bounding options.
    switch boundingMethod
        case 'dual-iter'
            options.nn.conzonotope_bounding_method = 'dual-iter';
            options.nn.conzonotope_bound_step_size = 1;
            options.nn.conzonotope_bound_max_iter = 1000;
        case 'fourier-motzkin'
            options.nn.conzonotope_bounding_method = 'fourier-motzkin';
            options.nn.polytope_bound_approx_max_iter = 8;
        case 'exact'
            options.nn.conzonotope_bounding_method = 'exact';
    end
    % Set default options.
    options = nnHelper.validateNNoptions(options);
    % Compute the bounds.
    [l,u,~,~] = conZonotope.approximateBoundsWithGPU(cZs,1,options);

    % [Optional] Visualize the bounds.
    % figure; hold on;
    % plot(zonotope(c,G),1:2,DisplayName='Zonotope');
    % plot(aux_2ConZonoWithEqConst(c,G,A,b),1:2,DisplayName='constrained Zonotope');
    % plot(interval(ldi,udi),1:2,DisplayName=sprintf('Bounds (%s)',boundingMethod));
    % legend
end

function cZ = aux_2ConZonoWithEqConst(c,G,A,b)
    % We convert the inequality constraints to equality constraints by 
    % adding a slack variable.

    % Obtain number of dimensions, generators, and batch size.
    [n,q] = size(G);
    % Obtain number of constraints.
    [p,~] = size(A);

    % Add a slack generators.
    G = [G zeros([n p])];
    % Compute scale for the slack variable.
    s = 1/2*(sum(abs(A),2) + b);
    A = [A eye(p).*s];
    % Compensate for the slack variable.
    b = b - s;
    
    % Instantiate constraint zonotope.
    cZ = conZonotope(c,G,A,b);
end

function cZs = aux_splitSet(nn,c,G,A,b)
    % Initialize cell array with current constraint zonotopes.
    cZs = {struct('c',c,'G',G,'A',A,'b',b)};
    % Iterate over the layers.
    for i=1:length(nn.layers)
        % Obtain i-th layer.
        layeri = nn.layers{i};
        % Check the type of layer.
        if isa(layeri,'nnLinearLayer')
            % Obtain the weight matrix and bias vector.
            Wi = layeri.W;
            bi = layeri.b;
            % Iterate all current sets.
            for j=1:length(cZs)
                % Obtain the j-th set.
                cZj = cZs{j};
                % Apply the affine map.
                cZs{j} = struct( ...
                    'c',Wi*cZj.c + bi,'G',Wi*cZj.G,'A',cZj.A,'b',cZj.b);
            end
        elseif isa(layeri,'nnReLULayer')
            % Obtain number of dimensions.
            nk = size(cZs{1}.c,1);
            % Iterate all dimensions.
            for k=1:nk
                % Obtain the number of current sets.
                numSplits = length(cZs);
                % Iterate all split sets.
                for j=1:numSplits
                    % Obtain the j-th set.
                    cZj = cZs{j};
                    % Compute slope of relu.
                    m = eye(nk);
                    m(k,:) = 0;
                    % <= 0
                    cZs{j} = struct('c',m*cZj.c,'G',m*cZj.G, ...
                        'A',[cZj.A; cZj.G(k,:)],'b',[cZj.b; -cZj.c(k)]);
                    % >= 0
                    cZs{end+1} = struct('c',cZj.c,'G',cZj.G, ...
                        'A',[cZj.A; -cZj.G(k,:)],'b',[cZj.b; cZj.c(k)]);
                end
            end
        else
            throw(CORAerror('CORA:nnLayerNotSupported', ...
                layeri,'splitInputSet'));
        end
    end
end

function f = aux_initFigure(Y)
    % Create the figure.
    f = figure; 
    % Create a subfigure for the input space.
    subplot(1,2,1); hold on; box on;
    title('Input Space')
    % Set the corrrect limits.
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])

    % Create a subfigure for the output space.
    subplot(1,2,2); hold on; box on;
    title('Output Space')
    % Set the correct limits.
    xlim([-1.75 -0.5])
    ylim([-3.25 -0.25])
    % Print the exact output set.
    plot(Y,1:2,'EdgeColor',CORAcolor('CORA:simulations'),'LineWidth',2)
end

function f = aux_plotSet(f,X,fid,varargin)
    % Plot a set X in the subfigure specified by fid.
    
    % Add a small interval to avoid numerical errors in plotting.
    pI = 1e-8*interval(-ones([ndims(X) 1]),ones([ndims(X) 1]));
    [color] = setDefaultValues({CORAcolor('CORA:reachSet')}, varargin);
    % Select the correct subplot.
    subplot(1,2,fid);
    % Plot the set.
    plot(X + pI,1:2,'FaceColor',color,'FaceAlpha',0.2,...
        'EdgeColor',color,'LineWidth',2);
end

% ------------------------------ END OF CODE ------------------------------
