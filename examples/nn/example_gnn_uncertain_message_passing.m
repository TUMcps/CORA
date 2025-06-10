function res = example_gnn_uncertain_message_passing
% example_gnn_uncertain_message_passing - example for enclosing a graph neural
%   networks with uncertain node features and uncertain graph structure
%
% Syntax:
%    res = example_gnn_uncertain_message_passing()
%
% Inputs:
%    -
%
% Outputs:
%    res - logical
%
% References:
%    [1] Ladner, T., et al. (2025). Formal Verification of Graph 
%        Convolutional Networks with Uncertain Node Features 
%        and Uncertain Graph Structure. TMLR.

% Authors:       Tobias Ladner
% Written:       17-February-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init initial set
X = interval( ...
    [0.9, 0.9; 1 1; 1 1], ...
    [1.1, 1.1; 1 1; 1 1] ...
);
X = reshape(X,[],1);
X = polyZonotope(X);

% init neural network
nn = neuralNetwork({ ...
   nnGCNLayer();  
   nnGCNLayer();
});
nn.reset();

% Evaluation with unknown presence of edge --------------------------------

% init graph
A = [
    1 1 1;
    1 1 1;
    1 1 1;
];
G = graph(A);

% network evaluation options
options = struct;
options.nn.graph = G;
options.nn.idx_pert_edges = 3; % Edge 1-3
options.nn.invsqrt_order = 2; % for larger graphs, 1 is sufficient

% evaluate 
timerVal = tic;
H_1 = nn.evaluate(X,options,1);
Y = nn.evaluate(H_1,options,2);
toc(timerVal);
drawnow

% obtain message passing
P = nnGCNLayer.message_passing.val;

% plot ---
disp("Plotting..")

figure; 

% plot graph
subplot(1,5,1); hold on; box on;
plot(G);
title('Graph')

% plot message passing
subplot(1,5,2); hold on; box on;
plot(P,2:3); 
plot(P.getSubset(7,1),2:3);
plot(P.getSubset(7,-1),2:3);
title('Message Passing');
xlabel('P_{(1,2)}')
ylabel('P_{(1,3)}')

% plot input
subplot(1,5,3); hold on; box on;
plot(X,[3,6]);
title('X')
xlabel('X_{(3)}')
ylabel('X_{(6)}')

% plot hidden layer
subplot(1,5,4); hold on; box on;
plot(H_1,[3,6]);
plot(H_1.getSubset(7,1),[3,6]);
plot(H_1.getSubset(7,-1),[3,6]);
title('H_1')
xlabel('H_{1(3)}')
ylabel('H_{1(6)}')

% plot output
subplot(1,5,5); hold on; box on;
plot(Y,[3,6], 'DisplayName','Enclosure uncertain (1)-(3)');
plot(Y.getSubset(7,1),[3,6],'DisplayName','Subset with (1)-(3)');
plot(Y.getSubset(7,-1),[3,6],'DisplayName','Subset without (1)-(3)');
title('Y')
xlabel('Y_{(3)}')
ylabel('Y_{(6)}')
legend


% Evaluation with edge ----------------------------------------------------

nn.resetGNN();

% network evaluation options
options = struct;
options.nn.graph = G;

% evaluate 
H_1 = nn.evaluate(X,options,1);
Y = nn.evaluate(H_1,options,2);

% obtain message passing
P = nnGCNLayer.message_passing.val;

% plot ---

subplot(1,5,2);
plot(P([2,2]),P([3,3]),'.');

subplot(1,5,4); hold on;
plot(H_1,[3,6]);

subplot(1,5,5); hold on;
plot(Y,[3,6], 'DisplayName','Exact with (1)-(3)');
legend


% Evaluation without edge -------------------------------------------------

nn.resetGNN();

% delete edge
A(1,3) = 0;
A(3,1) = 0;
G = graph(A);

% network evaluation options
options = struct;
options.nn.graph = G;

% evaluate 
H_1 = nn.evaluate(X,options,1);
Y = nn.evaluate(H_1,options,2);

% obtain message passing
P = nnGCNLayer.message_passing.val;

% plot ---

subplot(1,5,2);
plot(P([2,2]),P([3,3]),'.');

subplot(1,5,4); hold on;
plot(H_1,[3,6]);

subplot(1,5,5); hold on;
plot(Y,[3,6], 'DisplayName','Exact without (1)-(3)');
legend

% Compare P with interval matrices ----------------------------------------

A = interval(full(G.adjacency));
A(1,3) = interval(0,1);
A(3,1) = interval(0,1);
D_diag = sum(A);
P_int = diag(1./sqrt(D_diag)) * A * diag(1./sqrt(D_diag));

% plot
subplot(1,5,2);
plot(reshape(P_int,[],1),2:3,'--','DisplayName','Enclosure interval arithmetic')

% examples completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
