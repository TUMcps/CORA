function res = example_linear_verify_temporalLogic()
% example_linear_verify_temporalLogic - example for automated verification 
%       against temporal logic specifications for the HEAT benchmark from 
%       the ARCH competition using the algorithm from [1]
%
% Syntax:
%    res = example_linear_verify_temporalLogic()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% References:
%    [1] N. Kochdumper et al. "Fully-Automated Verification of Linear 
%        Time-Invariant Systems against Signal Temporal Logic 
%        Specifications via Reachability Analysis"

% Authors:       Niklas Kochdumper
% Written:       22-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------    

% load system matrices
load('heat10');

% construct output matric (center of block)
samples = nthroot(size(A,1),3);
indMid = 556;

C = zeros(1,size(A,1));
C(1,indMid) = 1;

% construct linear system object
sys = linearSys(A,B);


% Parameter ...............................................................
    
x0 = aux_getInit(A,samples,1);
temp = diag(0.1*x0);
temp = temp(:,x0 > 0);
params.R0 = zonotope(x0,temp);
params.U = zonotope(0);

params.tFinal = 40;

options = struct;
options.verifyAlg = 'stl:kochdumper';

% Specification -----------------------------------------------------------

x = stl('x',size(A,1));

eq = globally(x(151) > 0.07,interval(3,4)) | ...
                        finally(x(151) < 0.07,interval(4,6));
spec = specification(eq,'logic');


% Verification ------------------------------------------------------------

tic
[res,R,fals] = verify(sys, params, options, spec);
tComp = toc;

disp(['computation time for verification: ',num2str(tComp)]);
if res
    disp('specification satisfied!');
else
    disp('specification violated!');
end


% Visualization -----------------------------------------------------------

figure; hold on; box on;

% plot reachable set or falsifying trajectory
if res
    plotOverTime(R,151,'FaceColor',[.7 .7 .7]);
else
    simOpts.tFinal = params.tFinal;
    simOpts.x0 = fals.x0;

    sysDT = linearSysDT(sys,0.02);

    [t,x] = simulate(sysDT,simOpts);

    plot(t,x(:,151),'k');
end

% plot specifcation
plot([0,14],0.07*[1,1],'--r');
plot([3,3],[0,0.1],'--b');
plot([4,4],[0,0.1],'--b');

% formatting
xlim([0,6]); ylim([0,0.09]);
xlabel('time'); ylabel('x_{151}');
set(gcf, 'Units', 'centimeters', 'OuterPosition', [0, 0, 9, 8]);

end


% Auxiliary functions -----------------------------------------------------

function x0 = aux_getInit(A, samples, temp)
% returns an initial state vector

    x0 = zeros(size(A,1),1);

    if samples ~= 5 && (samples < 10 || mod(samples,10) ~= 0)
        throw(CORAerror('CORA:specialError',...
            'Init region is not evenly divided by discretization!'));
    end

    % maximum z point for initial region is 0.2 for 5 samples and 0.1 otherwise
    if samples >= 10
        max_z = (0.1/(1/samples));
    else
        max_z = (0.2/(1/samples));
    end

    for z = 0:max_z
        zoffset = z * samples * samples;

        for y = 0:(0.2/(1/samples))
            yoffset = y * samples;

            for x = 0:(0.4/(1/samples))
                index = x + yoffset + zoffset;

                x0(index+1) = temp;
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
