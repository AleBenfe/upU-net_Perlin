function [imdsTrain, imdsVal, imdsTest, pxdsTrain, pxdsVal, pxdsTest] = ...
    partitiondataIMAGES(imds,pxds,perc)
% Partition data by randomly selecting 60% of the data for training. The
% rest is used for testing.
    
% Set initial random state for example reproducibility.
% rng(0); 
numFiles = numel(imds.Files);
shuffledIndices = randperm(numFiles);

% Use 60% of the images for training.
N = round(perc(1) * numFiles);
% Use 50% of the images for validation and test
M = round((numFiles-N)/2)+1;

trainingIdx   = shuffledIndices(1:N);
validationIdx = shuffledIndices(N+1:(N+M));

% Use the rest for testing.
testIdx = shuffledIndices((N+M)+1:end);

% Create image datastores for training and test.
trainingImages = imds.Files(trainingIdx);
testImages = imds.Files(testIdx);
valImages = imds.Files(validationIdx);
imdsTrain = imageDatastore(trainingImages);
imdsTest = imageDatastore(testImages);
imdsVal = imageDatastore(valImages);


% Create pixel label datastores for training and test.
trainingLabels  = pxds.Files(trainingIdx);
testLabels      = pxds.Files(testIdx);
valLabels       = pxds.Files(validationIdx);
pxdsTrain       = imageDatastore(trainingLabels);
pxdsTest        = imageDatastore(testLabels);
pxdsVal         = imageDatastore(valLabels);
end