classdef levelSet < contSet
% levelSet - class definition for non-linear level sets
%
% Syntax:
%       obj = levelSet(eq,vars,compOp)
%       obj = levelSet(eq,vars,compOp,solved)
%
% Inputs:
%    eq - symbolic equation defining the set (eq == 0 or eq <= 0 or eq < 0)
%    vars - symbolic variables
%    compOp - operator ('==' or '<=' or '<')
%    solved - equation solved for one variable
%
% Outputs:
%    obj - generated level set object
%
% Example: 
%    syms x y
%    eq = x^2 + y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%    
%    figure
%    hold on
%    xlim([-3,3]);
%    ylim([-3,3]);
%    plot(ls,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: halfspace, constrainedHyperplane

% Author:       Niklas Kochdumper
% Written:      19-July-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    
    eq = [];        % symbolic equation
    compOp = [];    % comparison operator
    vars = [];      % symbolic variables
    funHan = [];    % function handle for the non-linear equation
    der = [];       % derivatives of the non-linear equation
    dim = [];       % dimension
    solved = [];    % symbolic equation solved for one variable
    solveable = 0;  % equation is solvable for one variable
    
end
    
methods
    
    % class constructor
    function obj = levelSet(eq,vars,compOp,varargin)
        
        if nargin == 1 && isa(eq,'levelSet')
            % copy constructor
            obj = eq;
            return
        end
        
        % set basic class properties
        obj.eq = eq;
        obj.vars = vars;
        obj.dim = length(vars);
        
        obj.funHan = matlabFunction(eq,'Vars',{vars});
        
        if ~iscell(compOp)
            if ismember(compOp,{'==','<=','<'})
               obj.compOp = compOp; 
            else
               error('levelSet: wrong value for input argument "compOp"!');
            end
        else
           for i = 1:length(compOp)
               if ~ismember(compOp,{'<=','<'})
                   error('levelSet: wrong value for input argument "compOp"!');
               end
           end
           obj.compOp = compOp;
        end
        
        % set parent object properties
        obj.dimension = obj.dim;
        
        % compute derivatives
        if strcmp(compOp,'==')
            grad = gradient(eq,vars);
            obj.der.grad =  matlabFunction(grad,'Vars',{vars});

            hess = hessian(eq,vars);
            obj.der.hess = matlabFunction(hess,'Vars',{vars});

            third = cell(length(grad),1);
            for i = 1:length(grad)
               temp = hessian(grad(i),vars); 
               third{i} = matlabFunction(temp,'Vars',{vars});
            end
            obj.der.third = third;
        end
        
        % try to solve non-linear equation for one variable
        if ~iscell(compOp) && strcmp(compOp,'==')
            
            if nargin > 3
                
                obj.solved = varargin{1};
                
            else
           
               obj.solved = cell(length(vars),1); 
               vars_ = symvar(obj.eq);

               % loop over all variables
               for i = 1:length(vars)

                  % check if variable is contained in non-linear equation
                  if ismember(vars(i),vars_)

                      obj.solved{i}.contained = 1;
                      temp = solve(eq == 0,vars(i));

                      % check if the equation could be solved for variable
                      try
                          if ~isempty(temp)
                             obj.solved{i}.solveable = 1;
                             obj.solved{i}.eq = temp;

                             % loop over all solutions
                             obj.solved{i}.funHan = cell(length(obj.solved{i}.eq),1);

                             for j = 1:length(obj.solved{i}.eq)

                                [eq_,grad_,hess_,third_] = derivatives(obj.solved{i}.eq(j),vars,i);

                                obj.solved{i}.funHan{j}.eq = eq_;
                                obj.solved{i}.funHan{j}.grad = grad_;
                                obj.solved{i}.funHan{j}.hess = hess_;
                                obj.solved{i}.funHan{j}.third = third_;

                                obj.solved{i}.funHan{j}.eq = matlabFunction(obj.solved{i}.eq(j),'Vars',{vars});
                             end

                             obj.solveable = 1;
                          else
                             obj.solved{i}.solveable = 0;
                             obj.solved{i}.eq = [];
                             obj.solved{i}.funHan = [];
                          end
                      catch
                          obj.solved{i}.solveable = 0;
                          obj.solved{i}.eq = [];
                      end
                  else
                      obj.solved{i}.contained = 0;
                      obj.solved{i}.solveable = 0;
                      obj.solved{i}.eq = [];
                      obj.solved{i}.funHan = [];
                  end
               end
            end
        end 
    end
end

end


% Auxiliary Functions -----------------------------------------------------

function [eq,grad,hess,third] = derivatives(eq,vars,i)
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

%------------- END OF CODE --------------