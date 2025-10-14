function [z] = computePGDAttack(obj,x,t,options,epsilon,Q,varargin)
% computePGDAttack - Compute a Projected Gradient Descent (PGD) adversarial attack.
%
% Syntax:
%    [z] = computePGDAttack(obj,x,y,options,epsilon,Q)
%
% Inputs:
%    obj - neural network
%    x - input value
%    t - target value
%    options - evaluation parameters (see neuralNetwork.evaluate)
%    epsilon - perturbation radius
%    Q - number of iterations
%    l,u - input bounds
%    stepsize - perturbation step size
%    stepsizeDecay - decay factor for stepsize
%    stepsizeDecayIter - iterations where decay factor is applied
%    lossDer - derivative of loss function
%    idxLayer - indices of layers to be evaluated
%
% Outputs:
%    z - altered input
%    
% References:
%    [1] A. Kurakin, I. J. Goodfellow, and S. Bengio. "Adversarial Machine 
%        Learning at Scale." In: 5th International Conference on Learning 
%        Representations ICLR. Toulon, France, 2017.
%    [2] A. Madry, A. Makelov, L. Schmidt, D. Tsipras, and A. Vladu. 
%        "Towards Deep Learning Models Resistant to Adversarial Attacks." 
%        In: 6th International Conference on Learning Representations 
%        (ICLR). Vancouver, BC, Canada, 2018.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: neuralNetwork

% Authors:       Lukas Koller
% Written:       21-June-2023
% Last update:   ---
% Last revision: ---    

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(6,13)
[l,u,stepsize,stepsizeDecay,stepsizeDecayIter,lossDer,idxLayer] = ...
    setDefaultValues({-inf,+inf,epsilon/Q,1,[],@(t,y) softmax(y) - t,...
        1:length(obj.layers)},varargin);

% validate input
inputArgsCheck({ ...
    {obj,'att','neuralNetwork'}; ...
    {x,'att','numeric'}; ...
    {t,'att','numeric'}; ... 
    {options,'att','struct'}; ... 
    {epsilon,'att','numeric'}; ...
    {Q,'att','numeric'};
})

l = max(l,x - epsilon);
u = min(u,x + epsilon);

z = x;
for q=1:Q
    if ismember(q,stepsizeDecayIter)
        stepsize = stepsize*stepsizeDecay;
    end
    % iterated FGSM attack
    z = computeFGSMAttack(obj,z,t,options,stepsize,lossDer,idxLayer);
    % clip values
    z = min(u,max(l,z));
end
end

% ------------------------------ END OF CODE ------------------------------
