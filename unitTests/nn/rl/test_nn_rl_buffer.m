function res = test_nn_rl_buffer()
% test_nn_rl_buffer - unit test function for
%     @buffer/buffer: instantiate buffer and add data. Test if maximum
%     length is exceeded.
%
% Syntax:
%    res = test_nn_rl_buffer()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Manuel Wendl
% Written:       27-August-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% give required options
options.rl.critic.nn.train.method = 'point';
options.rl.critic.nn.train.num_init_gens = 5;

b = buffer(10);
b = b.resetBuffer();

for i = 1:4
    b = b.fillBuffer({i*ones(2,1),i,i,i*ones(2,1),false},options);
end

% Test if buffer indx is correct
assert(b.currentIndx == 5);

b = b.resetBuffer();

% Test if buffer indx is reset correct
assert(b.currentIndx == 1);

for i = 1:11
    b = b.fillBuffer({rand(2,1),rand(1,1),rand(1,1),rand(2,1),false},options);
end

% Test if buffer indx is correct after exceeding maximal size
assert(b.currentIndx == 2);

% Test if buffer size is correct after more samples
assert(size(b.array{1},2) == b.bufferSize);

% Test set based
b = b.resetBuffer();
options.rl.critic.nn.train.method = 'set';
options.rl.critic.nn.train.num_init_gens = 5;

for i = 1:4
    b = b.fillBuffer({i*ones(2,1),zonotope(i,i),i,i*ones(2,1),false},options);
end

% Check dimension of stored zonotope
assert(size(b.array{2},2) == options.rl.critic.nn.train.num_init_gens + 1);

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
