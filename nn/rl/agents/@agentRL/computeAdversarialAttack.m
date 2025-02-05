function x_adv = computeAdversarialAttack(obj,x)
% computeAdversarialAttack - find adversarial attack
%
% Syntax:
%   x_adv = computeAdversarialAttack(obj,x)
%
% Inputs:
%   x - state of system
%
% Outputs:
%   x_adv - adversarial state 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: actor

% Authors:       Manuel Wendl
% Written:       05-December-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

x_adv = x;

% Simple Attacks
% random sampling from perturbation radius
if strcmp(obj.options.rl.actor.nn.train.method,'rand')
    factors = -1 + 2*rand(size(obj.options.rl.noise,2),1);
    x_adv = x + obj.options.rl.noise*factors;
end
% extreme sampling from perturbation radius
if strcmp(obj.options.rl.actor.nn.train.method,'extreme')
    factors = sign(-1 + 2*rand(size(obj.options.rl.noise,2),1));
    x_adv = x + obj.options.rl.noise*factors;
end

% Literature Attacks
% Naive attack:
% Attack is not implemented with a for loop as described in the paper but
% with a batch evaluation to speed up runtime. 
if strcmp(obj.options.rl.actor.nn.train.method,'naive')
    i = obj.options.rl.actor.nn.train.advOps.numSamples;
    alpha = obj.options.rl.actor.nn.train.advOps.alpha;
    beta = obj.options.rl.actor.nn.train.advOps.beta;

    epsilon = obj.options.rl.noise;

    a_star = obj.actor.nn.evaluate_(x,obj.options.rl.actor,obj.actor.idxLayer);
    Q_star = obj.targetCritic.nn.evaluate_(cat(1,x,a_star),obj.options.rl.critic,obj.targetCritic.idxLayer);

    n = betarnd(alpha,beta,[size(x,1),i]);

    xs_adv = x+n*epsilon;

    a_adv = obj.actor.nn.evaluate_(xs_adv,obj.options.rl.actor,obj.actor.idxLayer);
    Q_adv = obj.targetCritic.nn.evaluate(cat(1,xs_adv,a_adv),obj.options.rl.critic,obj.targetCritic.idxLayer);

    [minQ,ind] = min(Q_adv);
    if minQ < Q_star
        x_adv = xs_adv(:,ind);
    else
        x_adv = x;
    end
end

% Gradient based attack:
if strcmp(obj.options.rl.actor.nn.train.method,'grad')
    % read options
    obj.options.rl.critic.nn.train.updateGrad = false;

    i = obj.options.rl.actor.nn.train.advOps.numSamples;
    alpha = obj.options.rl.actor.nn.train.advOps.alpha;
    beta = obj.options.rl.actor.nn.train.advOps.beta;

    epsilon = obj.options.rl.noise;

    % evaluate
    a_star = obj.actor.nn.evaluate_(x,obj.options.rl.actor,obj.actor.idxLayer);
    Q_star = obj.targetCritic.nn.evaluate_(cat(1,x,a_star),obj.options.rl.critic,obj.targetCritic.idxLayer);
    % backprop
    grad = obj.targetCritic.nn.backprop(Q_star,obj.options.rl.critic,obj.targetCritic.idxLayer);
    grad = grad(1:size(x,1));
    grad_dir = grad/norm(grad);

    % compute adversary attack (TODO: FGSM?)
    n = betarnd(alpha,beta,[size(x,1),i]);
    xs_adv = x-n.*grad_dir*epsilon;

    % evaluate adversary example
    a_adv = obj.actor.nn.evaluate_(xs_adv,obj.options.rl.actor,obj.actor.idxLayer);
    Q_adv = obj.targetCritic.nn.evaluate(cat(1,xs_adv,a_adv),obj.options.rl.critic,obj.targetCritic.idxLayer);

    [minQ,ind] = min(Q_adv);
    if minQ < Q_star
        x_adv = xs_adv(:,ind);
    else
        x_adv = x;
    end

    obj.options.rl.critic.nn.train.updateGrad = false;
end

end

% ------------------------------ END OF CODE ------------------------------
