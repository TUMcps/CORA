close all; clear;

% equidistant time step vector
numPoints = 50;
timeStepVector = linspace(0,1,numPoints)';

% generate random data
data1 = 0.001*rand(1,numPoints);
data2 = 0.001*randn(1,numPoints);

% plot data
figure; hold on; box on;
h1 = plot(timeStepVector,data1);
h2 = plot(timeStepVector,data2);
xlabel('t');
ylabel('xi');
legend([h1,h2],'Data1','Data2','Location','northwest');

matlab2tikz('simpleplot.tikz');