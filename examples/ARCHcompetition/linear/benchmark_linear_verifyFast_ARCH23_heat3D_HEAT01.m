function text = benchmark_linear_verifyFast_ARCH23_heat3D_HEAT01
% benchmark_linear_verifyFast_ARCH23_heat3D_HEAT01 - heat3d benchmark from
%     the 2023 ARCH competition
%
% Syntax:
%    text = benchmark_linear_verifyFast_ARCH23_heat3D_HEAT01()
%
% Inputs:
%    -
%
% Outputs:
%    text - char array

% Authors:       Mark Wetzlinger
% Written:       23-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------    

% load system matrices
load('heat5');

% construct output matrix (center of block)
samples = nthroot(size(A,1),3);
indMid = ceil((size(A,1))/2);

C = zeros(1,size(A,1));
C(1,indMid) = 1;

% construct linear system object
sys = linearSys('heat',A,B,[],C);


% Parameters --------------------------------------------------------------
    
x0 = aux_getInit(A,samples,1);
temp = diag(0.1*x0);
temp = temp(:,x0 > 0);
params.R0 = zonotope(x0,temp);
params.U = zonotope(interval(zeros(size(sys.B,2),1)));

params.tFinal = 40;

options = struct();
options.verifyAlg = 'reachavoid:supportFunc';


% Specification -----------------------------------------------------------

% forall t: y1 <= 0.10379885
d = 0.10369885+1e-4;
spec = specification(halfspace(1,d),'safeSet');
% according to ARCH report, maximum 0.10369885 at time 9.44


% Verification ------------------------------------------------------------

% min steps needed: 1140
[res,fals,savedata] = verify(sys,params,options,spec);


% Return value ------------------------------------------------------------

text = ['Heat3D,HEAT01,',num2str(res),',',num2str(savedata.tComp)];

end


% Auxiliary functions -----------------------------------------------------

function x0 = aux_getInit(A, samples, temp)
% returns an initial state vector

    x0 = zeros(size(A,1),1);

    if samples ~= 5 && (samples < 10 || mod(samples,10) ~= 0)
       error('Init region is not evenly divided by discretization!')
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
