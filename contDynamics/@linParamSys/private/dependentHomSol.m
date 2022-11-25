function Rhom_tp = dependentHomSol(obj, Rinit, Uconst)
% dependentHomSol - computes the homogeneous solution when the parameters
% of the system matrix and the input are dependent
%
% Syntax:  
%    [obj] = dependentHomSol(obj, Rinit, Uconst)
%
% Inputs:
%    obj - linParamSys object 
%    Rinit - initial reachable set
%    Uconst - set of constant inputs
%
% Outputs:
%    Rhom_tp - homogeneous reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      02-June-2011
% Last update:  25-August-2011
% Last revision:---

%------------- BEGIN CODE --------------

%obtain required variables
Ac = obj.A.center;
Ag = obj.A.generator;
c = center(Rinit);
r = obj.stepSize;
dim = obj.dim;
Umat = Uconst.Z;
params = obj.A.gens;%'obj.A.gens' does not work anymore once lin_error2dAB adds lagrange remainder to system matrix

%SECOND ORDER DEPENDENT SOLUTION
%zero parametric order
R_c = c + Umat(:,1)*r;
for i=1:2
    R_c = R_c + Ac^i*r^i/factorial(i)*(c + Umat(:,1)*r/(i+1));
end

%first parametric order
%auxiliary value
M = eye(dim)*r + Ac*r^2/2 + Ac^2*r^3/6;
%loop
for i=1:params
    R_g(:,i) = (Ag{i}*r + Ac*Ag{i}*r^2/2 + Ag{i}*Ac*r^2/2) * c + ...
               (Ag{i}*r^2/2 + Ac*Ag{i}*r^3/6 + Ag{i}*Ac*r^3/6) * Umat(:,1) + ...
               M * Umat(:,i+1);
end

%second parametric order
%same index (i,i)
for i=1:params
    Rtmp = Ag{i}^2*r^2/2*c + ...
           Ag{i}*r^3/6*Umat(:,1) + ...
           (Ac*Ag{i} + Ag{i}*Ac)*r^3/6*Umat(:,i+1);
    R_g(:,end+1) = 0.5*Rtmp;
    R_c = R_c + 0.5*Rtmp;
end
%different index
if (params>=2)
    ind = combinator(params,2,'c');
    for i=1:length(ind(:,1))
        Atmp = Ag{ind(i,1)}*Ag{ind(i,2)} + Ag{ind(i,2)}*Ag{ind(i,1)};
        R_g(:,end+1) = Atmp*r^2/2*c + ...
                       Atmp*r^3/6*Umat(:,1) + ...
                       (Ac*Ag{ind(i,1)} + Ag{ind(i,1)}*Ac)*r^3/6*Umat(:,ind(i,2)+1);
    end
end
%obatin zonotope
R_lowOrder = zonotope([R_c,R_g]);

%HIGHER ORDER INDEPENDENT SOLUTION
%get state transition matrices
eZhigh = obj.mappingMatrixSet.highOrderZono;
eIhigh = obj.mappingMatrixSet.highOrderInt;
%get input transition matrices
eZhigh_input = obj.mappingMatrixSet.highOrderZonoInput;
eIhigh_input = obj.mappingMatrixSet.highOrderIntInput;

%remaining reachable set
R_rem_state = eZhigh*zonotope(c) + eIhigh*zonotope(c);
R_rem_input = eZhigh_input*Uconst + eIhigh_input*Uconst;
R_rem = R_rem_state + R_rem_input;

%STATE SOLUTION WITHOUT CENTER
Rinit_noCenter = Rinit + (-c);
R_hom_state = obj.mappingMatrixSet.zono*Rinit_noCenter + obj.mappingMatrixSet.int*Rinit_noCenter;

%FINAL SOLUTION
Rhom_tp = R_hom_state + R_lowOrder + R_rem;


%------------- END OF CODE --------------