function motivationPlot(plotInfo)
% motivationPlot - plots the motivation values for a lane change
%
% Syntax:
%    motivationPlot(plotInfo)
%
% Inputs:
%    ???
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       ???
% Written:       ???
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for i=1:length(plotInfo)
    %get vector for left lane motivation
    motLeft(i)=plotInfo{i}.left;
    %get vector for right lane motivation
    motRight(i)=plotInfo{i}.right;
    %get vector for following vehicle inconvenience
    motBehind(i)=plotInfo{i}.behind;
    %get vector for lane change probability
    prob(i)=plotInfo{i}.prob;    
end

hold on

%construct time vector
t=0.5:0.5:5;

%plot left lane motivation
plot(t,motLeft);
%plot right lane motivation
plot(t,motRight);
%plot following vehicle inconvenience
plot(t,motBehind);
%plot lane change probability
plot(t,prob);

% ------------------------------ END OF CODE ------------------------------
