classdef exponentialMatrix
% exponentialMatrix - class storing properties of the propagation matrix
%    e^At (incl. variations) as well as auxiliary values for re-usage
%
% Syntax:
%    obj = exponentialMatrix(A)
%
% Inputs:
%    A - system matrix (real-valued, square)
%
% Outputs:
%    obj - exponentialMatrix object
%
% Example:
%    A = [4 -1; 1 4];
%    expmat = exponentialMatrix(A);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       12-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
properties (SetAccess = private, GetAccess = public)
    A;                  % system matrix: A
    Aabs;               % absolute value of system matrix: |A|
    Apower;             % cell-array of A^i: A, A^2, A^3, ...
    Apower_abs;         % cell-array of |A|^i: |A|, |A|^2, |A|^3, ...

    Apos;               % cell-array of indices of positive values of A^i
    Aneg;               % cell-array of indices of negative values of A^i

    isAinv;             % true/false whether A is invertible
    Ainv;               % inverse of A: A^-1

    eADeltat;           % e^(A Delta t)
    eADeltat_back;      % e^(-A Delta t)

    eAtk;               % cell-array of matrices: e^At_1, e^At_2, ...
    eAtk_T;             % cell-array of transposed matrices: (e^At_1)^T, (e^At_2)^T, ...
    eAtk_back;          % cell-array of matrices: e^(-At_1), e^(-At_2), ...

    F;                  % curvature error
    Ftilde;             % curvature error (inputs)
    conv;               % true/false whether computation of F/Ftilde converged

end

methods
    function expmat = exponentialMatrix(A)

        % save system matrix
        expmat.A = A;
        

        % initialize first power of A and |A|
        expmat.Apower{1} = A;
        expmat.Apower_abs{1} = abs(A);
        % initialize positive and negative indices of A^eta
        expmat.Apos = cell(0);
        expmat.Aneg = cell(0);
        
        % pre-compute inverse of A matrix
        expmat.isAinv = rank(full(A)) == length(A);
        expmat.Ainv = [];
        if expmat.isAinv
            expmat.Ainv = inv(A);
        end
        

        % constructor call leaves some properties empty, use methods below:
        % - init_eADeltat
        expmat.eADeltat = [];
        expmat.eADeltat_back = [];

        % - init_eAtk, all_eAtk
        expmat.eAtk = cell(0);
        expmat.eAtk_T = cell(0);
        expmat.eAtk_back = cell(0);

        % - computeCurvatureErrors
        expmat.F = [];
        expmat.Ftilde = [];
        expmat.conv = [];

    end

    % class functions
    function expmat = init_eADeltat(expmat,dt)
    % initialize propagation matrices over a single step

        expmat.eADeltat = expm(expmat.A*dt);
        expmat.eADeltat_back = expm(-expmat.A*dt);
    end

    function expmat = init_eAtk(expmat,steps,varargin)
    % initializes the properties relating to the sequence of propagation
    % matrices e^Atk (including transpose and -A) to avoid costly
    % cell-array extension by +1 in each iteration

        % check if already initialized
        if ~isempty(expmat.eAtk) && ~isempty(expmat.eAtk_T) ...
                && ~isempty(expmat.eAtk_back)
            return
        end

        % start time
        tStart = setDefaultValues({0},varargin);

        % read out dimension
        n = length(expmat.A);

        % e^Atk
        expmat.eAtk = cell(steps+1,1);
        % transpose of e^Atk
        expmat.eAtk_T = cell(steps+1,1);
        % e^(-Atk)
        expmat.eAtk_back = cell(steps+1,1);

        % init first value
        if withinTol(tStart,0)
            expmat.eAtk{1} = eye(n);
            expmat.eAtk_T{1} = eye(n);
            expmat.eAtk_back{1} = eye(n);
        else
            expmat.eAtk{1} = expm(expmat.A*tStart);
            expmat.eAtk_T{1} = expmat.eAtk{1}';
            expmat.eAtk_back{1} = expm(-expmat.A*tStart);
        end
    end

    function expmat = next_Aposneg(expmat,eta)
    % the separation of A^eta into positive and negative indices can be
    % precomputed and saved; Apower_eta has to match eta correctly!
    
        if length(expmat.Apos) >= eta && ~isempty(expmat.Apos{eta})
            % ... then also length(expmat.Aneg) >= eta
            % already there (-> read from memory)
            
        else
            
            expmat.Aneg{eta} = expmat.Apower{eta};
            expmat.Apos{eta} = expmat.Apower{eta};
            expmat.Aneg{eta}(expmat.Aneg{eta} > 0) = 0;
            expmat.Apos{eta}(expmat.Apos{eta} < 0) = 0;
        end
    
    end

    function expmat = next_Apower(expmat,eta,varargin)
    % this function ensures that the eta-th power of A and |A| is computed
    % (this is necessary, since we do not know the largest power in
    % advance, and we want to save computations as much as possible)
    % we do not compute A^eta but A^eta / eta! instead to increase the
    % stability -> this has to be taken into account in all use cases!
    % (this is currently not enacted for |A|^eta)

    % choose whether Apower or Apower_abs should be computed
    flag = setDefaultValues({'standard'},varargin);
    inputArgsCheck({{expmat,'att','exponentialMatrix'}; ...
        {eta,'att','numeric',{'integer','positive'}}; ...
        {flag,'str',{'standard','abs'}}});
    
    
    switch flag
        case 'standard'
            % check A^eta
            if length(expmat.Apower) >= eta
                % read from memory
            else
                % compute all terms A^i/i! until eta
                maxeta = length(expmat.Apower);
                for i=maxeta:eta-1
                    expmat.Apower{i+1,1} = expmat.Apower{i} * expmat.A / (i+1);
                    % sparse/full  storage for more efficiency
%                     if nnz(expmat.Apower{i+1}) / (size(expmat.Apower{i+1},1)^2) < 0.1
%                         expmat.Apower{i+1} = sparse(expmat.Apower{i+1});
%                     else
                        expmat.Apower{i+1} = full(expmat.Apower{i+1});
%                     end
                end
            end

        case 'abs'
            % check |A|^eta
            if length(expmat.Apower_abs) >= eta
                % read from memory
            else
                % compute all powers |A|^i until eta
                maxeta = length(expmat.Apower_abs);
                for i=maxeta:eta-1
                    expmat.Apower_abs{i+1,1} = expmat.Apower_abs{i}*expmat.A_abs;
                    % sparse/full storage for more efficiency
%                     if nnz(expmat.Apower_abs{i+1}) / (size(expmat.Apower_abs{i+1},1)^2) < 0.1
%                         expmat.Apower_abs{i+1} = sparse(expmat.Apower_abs{i+1});
%                     else
                        expmat.Apower_abs{i+1} = full(expmat.Apower_abs{i+1});
%                     end
                end
            end

    end
    
    end

    function expmat = next_eAtk(expmat,k,varargin)
    % computes the next e^Atk, (e^Atk)^T or e^(-Atk)

        % choose which property should be computed
        flag = setDefaultValues({'standard'},varargin);
        inputArgsCheck({{expmat,'att','exponentialMatrix'}; ...
            {k,'att','numeric',{'integer','positive'}}; ...
            {flag,'str',{'standard','transpose','back'}}});

        switch flag
            case 'standard'
                % forward propagation
                if ~isempty(expmat.eAtk{k})
                    % already computed
                    return
                elseif isempty(expmat.eAtk{k-1})
                    % recursive call if (k-1)-th element not here
                    expmat = next_eAtk(expmat,k-1,flag);
                end
                expmat.eAtk{k} = expmat.eADeltat * expmat.eAtk{k-1};

            case 'transpose'
                % forward transpose (used for support function evaluation)
                if ~isempty(expmat.eAtk_T{k})
                    % already computed
                    return
                elseif isempty(expmat.eAtk{k})
                    % compute non-transposed variant
                    expmat = next_eAtk(expmat,k,'standard');
                end
                expmat.eAtk_T{k} = expmat.eAtk{k}';

            case 'back'
                % backward propagation
                if ~isempty(expmat.eAtk_back{k})
                    % already computed
                    return
                elseif isempty(expmat.eAtk_back{k-1})
                    % recursive call if (k-1)-th element not here
                    expmat = next_eAtk(expmat,k-1,flag);
                end
                expmat.eAtk_back{k} = expmat.eADeltat_back * expmat.eAtk_back{k-1};
        end

    end

    function expmat = computeCurvatureErrors(expmat,dt,varargin)
    % computation of F and Ftilde for a given time step size; currently,
    % these variables are computed by a Taylor series until floating-point
    % precision, i.e., we increase the truncation order until the
    % additional values are so small that the stored number (finite
    % precision!) does not change anymore

    % load data from object/options structure
    n = length(expmat.A);

    % check if input vector provided (otherwise assume 0)
    u = setDefaultValues({zeros(n,1)},varargin);
    
    % non-zero input vector? (skip computation of Ftilde if u is all-zero)
    isu = any(u);
    
    % initialize auxiliary variables and flags for loop
    Asum_pos_F = zeros(n);
    Asum_neg_F = zeros(n);
    stoploop_F = false;
    if isu
        Asum_pos_Ftilde = zeros(n);
        Asum_neg_Ftilde = zeros(n);
        stoploop_Ftilde = false;
    else
        expmat.Ftilde = [];
        stoploop_Ftilde = true;
    end
    
    eta = 1;
    while true
        % exponential: 1:eta
        
        % compute powers
        expmat = next_Apower(expmat,eta);
        
        % F starts at eta = 2, so skip for eta = 1
        if eta==1
            eta = eta + 1; continue
        end
        
        % tie/inputTie: 2:eta
        % note: usually, inputTie goes to eta+1 (with eta from F), but
        % since we compute terms until floating-point precision, this does
        % not need to be respected (only if we were to use a remainder term
        % E, which then would necessarily need to be adapted to a fixed eta)
        
        % compute factor (factorial already included in powers of A)
        exp1 = -(eta)/(eta-1);
        exp2 = -1/(eta-1);
        factor = ((eta)^exp1-(eta)^exp2) * dt^eta;
        
        if ~stoploop_F
            expmat = next_Aposneg(expmat,eta);
            
            % if new term does not change result anymore, loop to be finished
            Asum_add_pos_F = factor*expmat.Aneg{eta};
            Asum_add_neg_F = factor*expmat.Apos{eta};
            
            % safety check (if time step size too large, then the sum converges
            % too late so we already have Inf values)
            % previous criterion: any(any(isinf(Asum_add_pos_F))) || any(any(isinf(Asum_add_neg_F)))
            % ... but very costly for large A!
            if eta == 75
                expmat.F = [];
                expmat.Ftilde = [];
                expmat.conv = false;
                return
            end
            
            % compute ratio for floating-point precision
            if all(all(Asum_add_pos_F <= eps * Asum_pos_F)) && ...
                    all(all(Asum_add_neg_F >= eps * Asum_neg_F))
                stoploop_F = true;
            end
    
            % compute powers; factor is always negative
            Asum_pos_F = Asum_pos_F + Asum_add_pos_F; 
            Asum_neg_F = Asum_neg_F + Asum_add_neg_F;
    
        end
        
        if ~stoploop_Ftilde
            expmat = next_Aposneg(expmat,eta-1);
            
            % if new term does not change result anymore, loop to be finished
            % we require one additional division by eta as the terms in expmat
            % are divided by (eta-1)! instead of eta! as required
            Asum_add_pos_Ftilde = factor*expmat.Aneg{eta-1} / eta;
            Asum_add_neg_Ftilde = factor*expmat.Apos{eta-1} / eta;
            
            % safety check (if time step size too large, then the sum converges
            % too late so we already have Inf values)
            % previous criterion: any(any(isinf(Asum_add_pos_Ftilde))) || any(any(isinf(Asum_add_neg_Ftilde)))
            % ... but very costly for large A!
            if eta == 75
                expmat.F = [];
                expmat.Ftilde = [];
                expmat.conv = false;
                return
            end
            
            % compute ratio for floating-point precision
            if all(all(Asum_add_pos_Ftilde <= eps * Asum_pos_Ftilde)) && ...
                    all(all(Asum_add_neg_Ftilde >= eps * Asum_neg_Ftilde)) 
                stoploop_Ftilde = true;
            end
    
            % compute powers; factor is always negative
            Asum_pos_Ftilde = Asum_pos_Ftilde + Asum_add_pos_Ftilde; 
            Asum_neg_Ftilde = Asum_neg_Ftilde + Asum_add_neg_Ftilde;
        end
        
        % instantiate interval matrices if converged
        if stoploop_F
            expmat.F = interval(Asum_neg_F,Asum_pos_F);
        end
        if stoploop_Ftilde && isu
            expmat.Ftilde = interval(Asum_neg_Ftilde,Asum_pos_Ftilde);
        end
    
        % exit loop if both converged
        if stoploop_F && stoploop_Ftilde
            break;
        end
        
        % increment eta
        eta = eta + 1;
    end
    
    end
    
end

% methods (Static = true)
% 
% end

end

% ------------------------------ END OF CODE ------------------------------
