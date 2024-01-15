classdef levelSet < contSet
% levelSet - class definition for non-linear level sets
%
% Syntax:
%    obj = levelSet(eq,vars,compOp)
%    obj = levelSet(eq,vars,compOp,solved)
%
% Inputs:
%    eq - symbolic equation defining the set (eq == 0 or eq <= 0 or eq < 0)
%    vars - symbolic variables
%    compOp - operator ('==' or '<=' or '<')
%    solved - equation solved for one variable
%
% Outputs:
%    obj - generated levelSet object
%
% Example:
%    % single equation
%    syms x y
%    eqs = x^2 + y^2 - 4;
%    ls = levelSet(eqs,[x;y],'==');
%    
%    figure; hold on; xlim([-3,3]); ylim([-3,3]);
%    plot(ls,[1,2],'r');
%
%    % multiple equations
%    syms x y
%    eq1 = x^2 + y^2 - 4;
%    eq2 = x + y;
%    ls = levelSet([eq1;eq2],[x;y],{'<=';'<='});
%    
%    figure; hold on; xlim([-3,3]); ylim([-3,3]);
%    plot(ls,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace, conHyperplane

% Authors:       Niklas Kochdumper
% Written:       19-July-2019
% Last update:   ---
% Last revision: 16-June-2023 (MW, restructure using auxiliary functions)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    
    eq = [];            % symbolic equation
    vars = [];          % symbolic variables
    compOp = [];        % comparison operator
    solved = [];        % symbolic equation solved for one variable
    
    % internally-set properties
    funHan = [];        % function handle for the non-linear equation
    der = [];           % derivatives of the non-linear equation
    dim = [];           % dimension
    solvable = false;   % equation is solvable for one variable
    
end
    
methods
    
    % class constructor
    function obj = levelSet(varargin)

        % 0. avoid empty instantiation
        if nargin == 0
            throw(CORAerror('CORA:noInputInSetConstructor'));
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'levelSet')
            obj = varargin{1}; return
        end
        
        % 2. parse input arguments: varargin -> vars
        [eq,vars,compOp,solved] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(eq,vars,compOp,solved,nargin);

        % 4. compute internal properties
        [eq,vars,compOp,solved,funHan,der,dim,solvable] = ...
            aux_computeProperties(eq,vars,compOp,solved,nargin);

        % 5. assign properties
        obj.eq = eq;
        obj.vars = vars;
        obj.compOp = compOp;
        obj.solved = solved;
        obj.funHan = funHan;
        obj.der = der;
        obj.dim = dim;
        obj.solvable = solvable;
        
    end
end

methods (Static = true)
    ls = generateRandom(varargin) % generates a random level set
    ls = empty(n) % instantiates an empty level set
    ls = Inf(n) % instantiates a fullspace level set
end

end


% Auxiliary functions -----------------------------------------------------

function [eq,vars,compOp,solved] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 4
        throw(CORAerror('CORA:tooManyInputArgs',4));
    end

    % no input arguments
    if nargin == 0
        eq = []; vars = []; compOp = []; solved = [];
        return
    end
        
    % set basic class properties
    eq = varargin{1};
    vars = varargin{2};
    compOp = varargin{3};
    solved = [];
    if nargin == 4
        solved = varargin{4};
    end

end

function aux_checkInputArgs(eq,vars,compOp,solved,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        % check comparison operator
        if ~iscell(compOp)
            % single comparison operator has to be '==', '<=', '<'
            if ~ismember(compOp,{'==','<=','<'})
                throw(CORAerror('CORA:wrongValue','third',"be '==' or '<=' or '<'"));
            end
        else
            % for multiple comparison operators, no '==' allowed
            if any(cellfun(@(x) ~ismember(x,{'<=','<'}),compOp,'UniformOutput',true))
                throw(CORAerror('CORA:wrongValue','third',"be '<=' or '<'"));
            end
        end
    end

end

function [eq,vars,compOp,solved,funHan,der,dim,solvable] = aux_computeProperties(eq,vars,compOp,solved,n_in)
% compute properties according to given user inputs

    % set default values
    der = []; solvable = false; funHan = [];

    % dimension
    dim = length(vars);

    if n_in == 0
        return
    end
    
    % function handle
    funHan = matlabFunction(eq,'Vars',{vars});

    % compute derivatives
    if strcmp(compOp,'==')
        grad = gradient(eq,vars);
        der.grad =  matlabFunction(grad,'Vars',{vars});
    
        hess = hessian(eq,vars);
        der.hess = matlabFunction(hess,'Vars',{vars});
    
        third = cell(length(grad),1);
        for i = 1:length(grad)
            temp = hessian(grad(i),vars); 
            third{i} = matlabFunction(temp,'Vars',{vars});
        end
        der.third = third;
    end
    
    % try to solve non-linear equation for one variable
    if ~iscell(compOp) && strcmp(compOp,'==')
    
        if isempty(solved)
    
            solved = cell(length(vars),1); 
            vars_ = symvar(eq);
            
            % loop over all variables
            for i = 1:length(vars)
            
                % check if variable is contained in non-linear equation
                if ismember(vars(i),vars_)
                
                    solved{i}.contained = 1;
                    temp = solve(eq == 0,vars(i));
                
                    % check if the equation could be solved for variable
                    try
                        if ~isempty(temp)
                            solved{i}.solvable = 1;
                            solved{i}.eq = temp;
                            
                            % loop over all solutions
                            solved{i}.funHan = cell(length(solved{i}.eq),1);
                            
                            for j = 1:length(solved{i}.eq)
                                [eq_,grad_,hess_,third_] = aux_derivatives(solved{i}.eq(j),vars,i);
                                
                                solved{i}.funHan{j}.eq = eq_;
                                solved{i}.funHan{j}.grad = grad_;
                                solved{i}.funHan{j}.hess = hess_;
                                solved{i}.funHan{j}.third = third_;
                                
                                solved{i}.funHan{j}.eq = ...
                                matlabFunction(solved{i}.eq(j),'Vars',{vars});
                            end
                            
                            solvable = 1;

                        else
                            solved{i}.solvable = 0;
                            solved{i}.eq = [];
                            solved{i}.funHan = [];
                        end
                    catch
                        solved{i}.solvable = 0;
                        solved{i}.eq = [];
                    end

                else
                    solved{i}.contained = 0;
                    solved{i}.solvable = 0;
                    solved{i}.eq = [];
                    solved{i}.funHan = [];
                end
            end
        end
    end 

end

function [eq,grad,hess,third] = aux_derivatives(eq,vars,i)
% comptue derivatives and convert them to function handles
    
    % remove current variable
    vars_ = vars;
    vars_(i) = [];

    % equation
    eq_ = eq;
    eq = matlabFunction(eq_,'Vars',{vars});
    
    % gradient
    grad_ = gradient(eq_,vars_);
    grad =  matlabFunction(grad_,'Vars',{vars});
    
    % hessian matrix
    hess_ = hessian(eq_,vars_);
    hess = matlabFunction(hess_,'Vars',{vars});
    
    % third-order tensor
    third = cell(length(grad_),1);
    for i = 1:length(grad_)
       temp = hessian(grad_(i),vars_); 
       third{i} = matlabFunction(temp,'Vars',{vars});
    end
end

% ------------------------------ END OF CODE ------------------------------
