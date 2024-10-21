function res = testLong_pointSegmentSignal_set
% testLong_pointSegmentSignal_set - long unit test function of set
%
% Syntax:
%    res = testLong_pointSegmentSignal_set
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Florian Lercher
% Written:       16-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parameters
nrTests = 1000;
nrSamples = 10;
randFactor = 100;

% value patterns
tt = true;
ff = false;
pattern(1,:) = [tt,tt,ff,ff,tt,tt,ff,ff,tt,tt];
pattern(2,:) = [ff,ff,tt,tt,ff,ff,tt,tt,ff,ff];
pattern(3,:) = [tt,ff,tt,ff,tt,ff,tt,ff,tt,ff];
pattern(4,:) = [ff,tt,ff,tt,ff,tt,ff,tt,ff,tt];
pattern(5,:) = [tt,tt,tt,ff,tt,ff,tt,ff,tt,ff];
nrTimePoints = size(pattern,2) / 2;
nrSignals = size(pattern,1);

for i = 1:nrTests
    % generate random time points
    tp(1) = 0;
    tp(2:nrTimePoints) = sort(rand(1,nrTimePoints - 1)) .* randFactor;

    % generate random interval bounds
    bounds = sort(rand(1,2)) .* randFactor;
    
    % for each value pattern ...
    for j = 1:nrSignals
        sig = pointSegmentSignal(tp,pattern(j,:));
        % sample random points
        randPoints = rand(1,nrSamples) .* randFactor;
        samples = [randPoints,bounds];
        % for each open/closed combination ... 
        for k = 0:3
            lc = mod(floor(k/2),2) == 0;
            rc = mod(k,2) == 0;
            int = stlInterval(bounds(1),bounds(2),lc,rc);
            sigT = sig.set(int,tt);
            sigF = sig.set(int,ff);
            % for each sample ...
            for l = 1:length(samples)
                point = samples(l);
                % check if the signal is correct
                if contains(int,point)
                    assertLoop(sigT.at(point) || sigF.at(point),i,j)
                else
                    assertLoop(sigT.at(point) == sig.at(point),i,j)
                    assertLoop(sigF.at(point) == sig.at(point),i,j)
                end
            end
        end
    end
end

res = true;
end

% ------------------------------ END OF CODE ------------------------------
