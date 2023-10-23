function [obj,options] = expmtie_adaptive(obj,options)
% expmtie_adaptive - computes the remainder of the exponential matrix
%    and the correction matrix
%
% Syntax:
%    obj = expmtie_adaptive(obj,options)
%
% Inputs:
%    obj - linearSys object
%    options - options struct
%
% Outputs:
%    obj - linearSys object
%    options - options struct
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% References: 
%   -
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       26-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load data from object/options structure
A = obj.A;
A_abs = abs(A);
n = obj.dim;
deltat = options.timeStep;
% deltatbyfac = options.factor;

% initialize 
Apower{1} = A;
Apower_abs{1} = A_abs;
M = eye(n);
Asum_pos = zeros(n);
Asum_neg = zeros(n);
    
% increment taylorTerms until convergence
eta = 1;
while true
    % exponential: 1:eta
    % compute powers
    Apower{eta+1} = Apower{eta}*A;
    Apower_abs{eta+1} = Apower_abs{eta}*A_abs;
    deltatbyfac(eta,1) = deltat^eta / factorial(eta);
    M = M + Apower_abs{eta}*deltatbyfac(eta);
    
    if eta==1; eta = eta+1; continue; end % F starts at eta=2
    
    % tie: 2:eta
    % compute factor
    exp1 = -(eta)/(eta-1); exp2 = -1/(eta-1);
    factor = ((eta)^exp1-(eta)^exp2) * deltatbyfac(eta); 
    
    % init Apos, Aneg
    Apos = zeros(n);
    Aneg = zeros(n);
    
    % obtain positive and negative parts
    pos_ind = Apower{eta} > 0;
    neg_ind = Apower{eta} < 0;
    
    Apos(pos_ind) = Apower{eta}(pos_ind);
    Aneg(neg_ind) = Apower{eta}(neg_ind);
    
    % compute powers; factor is always negative
    Asum_pos = Asum_pos + factor*Aneg; 
    Asum_neg = Asum_neg + factor*Apos;
    
    % norm as stopping criterion (no abs() because already all positive)
    normAfro(eta,1) = norm(Asum_pos - Asum_neg,'fro');
    
    if ~any(any(Apower{eta})) ... % nilpotent
            || 1 - normAfro(eta-1)/normAfro(eta) < options.zetaTlin % relative convergence
        % determine error due to finite Taylor series, see eq.(2) in [1]
        W = expm(A_abs*options.timeStep) - M;
        % compute absolute value of W for numerical stability
        W = abs(W);
        E = interval(-W,W);
        % instantiate interval matrix
        Asum = interval(Asum_neg,Asum_pos);
        deltatbyfac(eta+1,1) = deltat^(eta+1) / factorial(eta+1); %inputSolution
        options.taylorTerms = eta;
        options.factor = deltatbyfac;
        options.etalinFro = 1 - normAfro(eta-1)/normAfro(eta);
        break
    else
        eta = eta + 1;
        if eta > 50
            throw(MException('expmtie:notconverging',...
                'Time Step Size too big.'));
        end
    end
end

% write to object structure
obj.taylor.powers = Apower;
obj.taylor.error = E;
obj.taylor.F = Asum + E;

end

% ------------------------------ END OF CODE ------------------------------
