function randBatch = getRandomBatch(obj,options)
% getRandomBatch - get a randomly sampled batch from the replay buffer
%
% Syntax:
%   randBatch = getRandomBatch(obj,options)
%
% Inputs:
%   options.rl - reinforcement learning options
% 
% Outputs:
%   randBatch - random training batch from replpay buffer
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: buffer

% Authors:       Manuel Wendl
% Written:       03-November-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

inputArgsCheck({ ...
    {options, 'att', 'struct'}
    })
if options.rl.batchsize > obj.currentIndx
    throw(CORAerror('CORA:specialError','Demanded batchsize exceeds the filling of the replay buffer.'))
end

ind = randperm(obj.currentIndx-1,options.rl.batchsize);
randBatch = cell(1,5);

randBatch{1} = obj.array{1}(:,ind);

if strcmp(options.rl.critic.nn.train.method,'set')
    randBatch{2} = obj.array{2}(:,:,ind);
else
    randBatch{2} = obj.array{2}(:,ind);
end
randBatch{3} = obj.array{3}(:,ind);
randBatch{4} = obj.array{4}(:,ind);
randBatch{5} = obj.array{5}(:,ind);
end

% ------------------------------ END OF CODE ------------------------------
