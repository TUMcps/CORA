function res = example_neuralNetwork_relu_inputTiling()
% example_neuralNetwork_relu_inputTiling - example for the visualization of 
%   the polytope-tiling of the input space by a ReLU-neural network.
%
% Syntax:
%    res = example_neuralNetwork_relu_inputTiling()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% References: -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Lukas Koller
% Written:       09-July-2025
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

% Obtain number of layers.
kappa = length(nn.layers);

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
C = [-2 1]; 
d = [-2.6];
U = polytope(C,d);

% Initialize input set.
Xi = X;

% Initialize figure;
f = aux_initFigure(Y);
% Plot initial input set.
f = aux_plotInput(f,Xi);
% Plot specification.
f = aux_plotOutput(f,U,CORAcolor('CORA:unsafe'));

% Sample input points (along each dimension).
N = 1e6;
xs = randPoint(zonotope(cx,Gx),N,'uniform');

% Compute the activation pattern for each input.
options.nn.train.backprop = true;
ys = nn.evaluate(xs,options);

% Iterate over the layers an extract the activation pattern.
actPats = [];

for i=1:length(nn.layers)
    % Obtain the i-th layer.
    layeri = nn.layers{i};
    if isa(layeri,'nnReLULayer')
        % Obtain the stored input.
        xi = layeri.backprop.store.input;
        % Extract the activation pattern.
        actPati = (xi > 0);
        % Append the activation pattern.
        actPats = [actPats; actPati];
    end
end

% Initialize cell for storing the polytopes.
Ps = {};

% Initialize an array for points that have not been added to a polytope
% yet.
hasBeedAdded = zeros(1,N,'logical');

% Construct the polytope for each activation pattern.
while ~all(hasBeedAdded)
    % Obtain the first activation pattern that has not beed added yet.
    idx = find(~hasBeedAdded,1);
    actPati = actPats(:,idx);
    % Find all points with the same activation pattern.
    ids = all(actPats == actPati,1);
    % Construct the polytope.
    xsi = xs(:,ids);
    k = convhull(xsi');
    % Append the polytope.
    Ps = [Ps {polytope(xsi(:,k))}];
    % Mark all samples as visited.
    hasBeedAdded(ids) = true;
end

% Plot all polytopes.
useCORAcolors('CORA:manual');
for i=1:length(Ps)
    % Obtain the color.
    pltClri = CORAcolor('CORA:next');
    % Plot the polytope in the input space.
    f = aux_plotInput(f,Ps{i},pltClri);
    % Pick any point from the polytope to compute the affine map.
    xi = randPoint(Ps{i},1);
    % Initialize the affine transform.
    Wi = diag(ones(size(xi)));
    bi = zeros(size(xi));
    for j=1:length(nn.layers)
        % Obtain the j-th layer.
        layerj = nn.layers{j};
        if isa(layerj,'nnReLULayer')
            % Obtain the activation pattern.
            actPatj = (xi > 0);
            % Compute the output.
            xi = max(xi,0);
            % Compute the resulting affine transform.
            Wi = diag(actPatj)*Wi;
            bi = diag(actPatj)*bi;
        elseif isa(layerj,'nnLinearLayer')
            % Obtain the affine transform of the j-th layer.
            Wj = layerj.W;
            bj = layerj.b;
            % Compute the output.
            xi = Wj*xi + bj;
            % Compute the resulting affine transform.
            Wi = Wj*Wi;
            bi = Wj*bi + bj;
        else
            throw(CORAerror('CORA:notSupported',...
                sprintf("Layer not supported '%s'!",layerj.name)));
        end
    end
    % Compute the output set for the input polytope by applying the affine
    % transform.
    Yi = Wi*Ps{i} + bi;
    % Plot the output polytope in the output space.
    f = aux_plotOutput(f,Yi,pltClri);
end

end


% Auxiliary functions -----------------------------------------------------

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
            % Current only linear and relu layers are supported.
            throw(CORAerror('CORA:nnLayerNotSupported', ...
                layeri,'splitInputSet'));
        end
    end
end

function f = aux_initFigure(Y)
    % Initialize a figure.
    f = figure; 
    % Create the input-subplot.
    subplot(1,2,1); hold on; box on;
    % Set title.
    title('Input Space')
    % Set axis limits.
    xlim([-1.25 1.25])
    ylim([-1.25 1.25])
    % Create output subplot.
    subplot(1,2,2); hold on; box on;
    % Set title.
    title('Output Space')
    % Specify axis limits.
    xlim([-1.75 -0.5])
    ylim([-3.25 -0.25])
    % Plot the output set.
    plot(Y,1:2,'EdgeColor',CORAcolor('CORA:simulations'),'LineWidth',2)
end

function f = aux_plotInput(f,X,varargin)
    % Create an interval to avoid numerical issues.
    pI = 1e-8*interval(-ones([ndims(X) 1]),ones([ndims(X) 1]));
    % Parse arguments.
    [color,faceAlpha] = setDefaultValues( ...
        {CORAcolor('CORA:reachSet'),0.2},varargin);
    % Select the correct subplot for the input space.
    subplot(1,2,1);
    % Plot the given set in the input space.
    plot(X + pI,1:2,'FaceColor',color,'FaceAlpha',faceAlpha,...
        'EdgeColor',color,'LineWidth',2);
end

function f = aux_plotOutput(f,Y,varargin)
    % Create an interval to avoid numerical issues.
    pI = 1e-8*interval(-ones([ndims(Y) 1]),ones([ndims(Y) 1]));
    % Parse arguments.
    [color,faceAlpha] = setDefaultValues( ...
        {CORAcolor('CORA:reachSet'),0.2},varargin);
    % Select the correct subplot for the output space.
    subplot(1,2,2);
    % Plot the given set in the output space.
    plot(Y + pI,1:2,'FaceColor',color,'FaceAlpha',faceAlpha, ...
        'EdgeColor',color,'LineWidth',2);
end


% ------------------------------ END OF CODE ------------------------------
