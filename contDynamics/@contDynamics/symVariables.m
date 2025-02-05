function [vars,vars_der] = symVariables(sys,varargin)
% symVariables - generates symbolic variables of a continuous system; the
%    symbolic variables are either set or a 0x1 sym object
%
% Syntax:
%    [vars,vars_der] = symVariables(sys)
%    [vars,vars_der] = symVariables(sys,withBrackets)
%
% Inputs:
%    sys - contDynamics object
%    withBrackets - true/false
%       true:  variable 'x' results in symbolic variables x1, x2, ...
%       false: variable 'x' results in symbolic variables xL1R, xL2R, ...
%
% Outputs:
%    vars - struct with fields
%       .x - symbolic state variables
%       .u - symbolic input variables
%       .y - symbolic constraint variables
%       .o - symbolic output variables
%       .p - symbolic parameters
%    vars_der - struct with fields
%       .dx - symbolic state deviation from linearization point
%       .du - symbolic input deviation from linearization point
%       .dy - symbolic constraint deviation from linearization point
%       .do - symbolic output deviation from linearization point
%
% Example: 
%    sys = contDynamics('test',3,1,2);
%    [vars,vars_der] = symVariables(sys);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: derivatives

% Authors:       Matthias Althoff, Mark Wetzlinger, Tobias Ladner
% Written:       18-January-2008
% Last update:   06-July-2017
%                05-November-2017
%                14-January-2018
%                19-November-2022 (MW, add outputs)
%                03-February-2023 (SM, real symbolic variables)
%                05-March-2024 (LL, consider nonlinearARX sys)
%                07-October-2024 (MW, use auxiliary function)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(1,2);
withBrackets = setDefaultValues({false},varargin);

%%% TODO: fix possibly wrong usage of sys.nrOfDims for nonlinearARX and
%%%       harmonize calls below

% symbolic variables
vars = struct('x',[],'u',[],'y',[],'o',[],'p',[]);
% symbolic variables for deviation from linearization point
vars_der = struct('x',[],'u',[],'y',[],'o',[]);

% states and inputs
if isa(sys,'nonlinearARX')
    vars.x = aux_symVector('x',sys.n_p*sys.nrOfOutputs,withBrackets);
    vars_der.x = aux_symVector('dx',sys.n_p*sys.nrOfOutputs,withBrackets);
    vars.u = aux_symVector('u',(sys.n_p+1)*sys.nrOfInputs,withBrackets);
    vars_der.u = aux_symVector('du',(sys.n_p+1)*sys.nrOfInputs,withBrackets);
else
    vars.x = aux_symVector('x',sys.nrOfDims,withBrackets);
    vars_der.x = aux_symVector('dx',sys.nrOfDims,withBrackets);
    vars.u = aux_symVector('u',sys.nrOfInputs,withBrackets);
    vars_der.u = aux_symVector('du',sys.nrOfInputs,withBrackets);
end

% algebraic constraints
if isprop(sys,'nrOfConstraints')
    vars.y = aux_symVector('y',sys.nrOfConstraints,withBrackets);
    vars_der.y = aux_symVector('dy',sys.nrOfConstraints,withBrackets);
else
    vars.y = sym('y',[0,1]);
    vars_der.y = sym('dy',[0,1]);
end

% outputs
vars.o = aux_symVector('y',sys.nrOfOutputs,withBrackets);
vars_der.o = aux_symVector('do',sys.nrOfOutputs,withBrackets);

% parameters
if isprop(sys,'nrOfParam')
    vars.p = aux_symVector('p',sys.nrOfParam,withBrackets);
else
    vars.p = sym('p',[0,1]);
end

end


% Auxiliary functions -----------------------------------------------------

function symVars = aux_symVector(varName,nrVars,withBrackets)
% withBrackets == true sandwiches the number in between 'L' and 'R'

if withBrackets
    % apparently, this case differentiation below is necessary...
    if nrVars == 1
        symVars = sym([varName 'L%dR'],[2,1],'real');
        symVars = symVars(1);
    else
        symVars = sym([varName 'L%dR'],[nrVars,1],'real');
    end
else
    symVars = sym(varName,[nrVars,1],'real');
end

end


% ------------------------------ END OF CODE ------------------------------
