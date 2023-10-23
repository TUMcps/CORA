function res = testLong_linParamSys_F
% testLong_linParamSys_F - underapproximate expm(A) and F for an interval matrix using
% sampling; check if samples are contained in over-approximation.
%
% Syntax:
%    res = testLong_linParamSys_F
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false if successfully completed
%

% Authors:       Matthias Althoff
% Written:       19-January-2007
% Last update:   28-July-2020
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enable access to private function "tie" and "mappingMatrix"
path = CORAROOT;
source1 = fullfile(path,'contDynamics','@linParamSys','private','tie.m');
target1 = fullfile(path,'contDynamics','@linParamSys','tie.m');
copyfile(source1,target1);
source2 = fullfile(path,'contDynamics','@linParamSys','private','mappingMatrix.m');
target2 = fullfile(path,'contDynamics','@linParamSys','mappingMatrix.m');
copyfile(source2,target2);
rmpath(genpath(path));
addpath(genpath(path));

%instantiate linear dynamics with constant parameters
S = rand(2); % random single value matrix
deltaS = rand(2); % random delta matrix
Aint = intervalMatrix(S, deltaS);
linSys  = linParamSys(Aint, 1);

% compute overapproximative F matrix
options.intermediateTerms = 2;
options.timeStep = 1;
mappingMatrix(linSys,options); % required for tie method
linSys = tie(linSys);

% extract interval matrix of exponential matrices
expm_int = linSys.mappingMatrixSet.zono + linSys.mappingMatrixSet.int;

% choose step sizes for time and system matrix
times = 10;
r = options.timeStep;
delta_t = r/times;
segments = [4,4,4,4];
nrOfSegments = prod(segments);
dim = length(S);

%init Exp_min, Exp_max, A
Expm_min=expm(S*r);
Expm_max=expm(S*r);
A=S;
A_min=S; A_max=S;


% get stepS (change of S in each iteration)
stepS = S; % preallocation
for iPos=1:length(segments)
    stepS(iPos)=2*deltaS(iPos)/(segments(iPos)-1);
end

% iterate over system matrices A
for i=1:nrOfSegments

    seg=i2s(segments,i); 
    for iPos=1:length(seg)
        A(iPos) = S(iPos) - deltaS(iPos) + (seg(iPos)-1)*stepS(iPos);
    end

    % compute exponential matrix for final time
    Exp_r = expm(A*r);
        
    % iterate over times
    for iTime=1:times
        t=delta_t*iTime;
    
        % compute F
        Exp_t = expm(A*t);
        F = Exp_t - eye(dim) - t/r*(Exp_r-eye(dim));
        
        % find minimum and maximum F and A
        if iTime==1 && i==1
            F_min=F;
            F_max=F;
        else
            for iPos=1:length(seg)
                
                %F_min/max
                if F(iPos)<F_min(iPos)
                    F_min(iPos) = F(iPos);
                elseif F(iPos)>F_max(iPos)
                    F_max(iPos) = F(iPos);
                end
                
                %Exp_min/max
                if Exp_r(iPos)<Expm_min(iPos)
                    Expm_min(iPos) = Exp_r(iPos);
                elseif Exp_r(iPos)>Expm_max(iPos)
                    Expm_max(iPos) = Exp_r(iPos);
                end
                
                %A_min/max
                if A(iPos)<A_min(iPos)
                    A_min(iPos) = A(iPos);
                elseif A(iPos)>A_max(iPos)
                    A_max(iPos) = A(iPos);
                end
            end
        end
    end
end

% delete copied files
delete(target1);
delete(target2);
rmpath(genpath(path));
addpath(genpath(path));

% check if lower bound of F included
res(1) = all(all(infimum(linSys.F.int) <= F_min));

% check if upper bound of F included
res(2) = all(all(supremum(linSys.F.int) >= F_max));

% check if lower bound of Expm included
res(3) = all(all(infimum(expm_int.int) <= Expm_min));

% check if upper bound of Expm included
res(4) = all(all(supremum(expm_int.int) >= Expm_max));

% check if lower and upper bound holds
res = all(res);

% ------------------------------ END OF CODE ------------------------------
