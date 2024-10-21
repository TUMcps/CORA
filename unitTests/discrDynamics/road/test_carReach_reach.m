function res = test_carReach_reach
% test_carReach_reach - unit test function for testing probabilistic
%    prediction of road vehicles using reachability analysis
%
%    Checks the generation of a small Markov chain for probabilistic
%    prediction of traffic participants; The Markov model is generated
%    based on reachability analysis
%
% Syntax:
%    res = test_carReach_reach
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Matthias Althoff
% Written:       31-July-2016
% Last update:   31-July-2017
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set fileName and pathName
fileName = 'fArray_unitTest.mat';
pathName = [CORAROOT filesep 'discrDynamics' filesep 'ProbOccupancyPrediction'];
modelInitialization = @initCar_unitTest;

% compute probabilistic model
probModel = carReach_reach(fileName,pathName,modelInitialization);

% obtain results-----------------------------------------------------------
T = probModel.T;
projMat = probModel.projMat;
GammaFull = probModel.GammaFull;
finalTime = probModel.timeStep;
stateField = probModel.stateField;
inputField = probModel.inputField;
%--------------------------------------------------------------------------

% load ground-truth results------------------------------------------------
load probModel_reach_groundTruth

T_groundTruth = probModel_groundTruth.T;
projMat_groundTruth = probModel_groundTruth.projMat;
GammaFull_groundTruth = probModel_groundTruth.GammaFull;
finalTime_groundTruth = probModel_groundTruth.timeStep;
stateField_groundTruth = probModel_groundTruth.stateField;
inputField_groundTruth = probModel_groundTruth.inputField;
%--------------------------------------------------------------------------

% perform comparisons------------------------------------------------------
assert(max(max(abs(T.T - T_groundTruth.T))) < 1e-12);
assert(max(max(abs(T.OT - T_groundTruth.OT))) < 1e-12);
assert(max(max(abs(projMat - projMat_groundTruth))) < 1e-12);
assert(max(max(abs(GammaFull - GammaFull_groundTruth))) < 1e-12);
assert(abs(finalTime - finalTime_groundTruth) < 1e-12);
%state field
assert(max(max(abs(stateField.intervals - stateField_groundTruth.intervals))) < 1e-12);
assert(max(abs(stateField.nrOfSegments - stateField_groundTruth.nrOfSegments)) < 1e-12);
segmentLength = (stateField.intervals(:,2)-stateField.intervals(:,1))./stateField.nrOfSegments;  % <-- segment Length is no longer a field AP
groundTruthSegmentLength = (stateField_groundTruth.intervals(:,2)-stateField_groundTruth.intervals(:,1))./stateField_groundTruth.nrOfSegments;  % <-- segment Length is no longer a field AP
assert(max(abs(segmentLength - groundTruthSegmentLength)) < 1e-12);
%input field
assert(max(max(abs(inputField.intervals - inputField_groundTruth.intervals))) < 1e-12);
assert(max(abs(inputField.nrOfSegments - inputField_groundTruth.nrOfSegments)) < 1e-12);
segmentLength = (inputField.intervals(:,2)-inputField.intervals(:,1))./inputField.nrOfSegments;  % <-- segment Length is no longer a field AP
groundTruthSegmentLength = (inputField_groundTruth.intervals(:,2)-inputField_groundTruth.intervals(:,1))./inputField_groundTruth.nrOfSegments;  % <-- segment Length is no longer a field AP
assert(max(abs(segmentLength - groundTruthSegmentLength)) < 1e-12);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
