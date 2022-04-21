clearvars
close all
clc

%Select the problem size 
imgSize = 128;

% Select if test a pretrained network (0) or train a new one (>0)
training = 0;

% Create the data stores 
fullSet = imageDatastore(sprintf('Images/%d',imgSize));
fullRes = imageDatastore(sprintf('Labels/%d',imgSize));

% Divide the dataset in training, validation and test sets. 
[imdsTrain, imdsTest, imdsVal,...
    pxdsTrain, pxdsTest, pxdsVal] = partitiondataIMAGES(fullSet,fullRes,0.9);


dsTrain = combine(imdsTrain,pxdsTrain);
dsVal   = combine(imdsVal,pxdsVal);
dsTest  = combine(imdsTest,pxdsTest);


if training
    % Create the Network structure
    inputSize = [size(imread(imdsTrain.Files{1})), 1];
    convlayers=setNetwork(inputSize,4,8);

    % Set the training options
    options = trainingOptions('adam', ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropFactor',0.2, ...
        'LearnRateDropPeriod',25, ...
        'L2Regularization',0.1,...
        'MaxEpochs',200, ...
        'MiniBatchSize',50, ...
        'InitialLearnRate',1e-1, ...
        'Shuffle','every-epoch', ...
        'Plots','training-progress', ...
        'ExecutionEnvironment','parallel',...
        'ValidationData',dsVal,...
        'ValidationFrequency',150,...
        'Verbose',true);

    % Train the network
    net = trainNetwork(dsTrain,convlayers,options);
else
    % Load previously trained network 
    if size(imread(imdsTrain.Files{1}),1)==128
        load upUNet_128.mat
    elseif size(imread(imdsTrain.Files{1}),1)==256
        load upUNet_256.mat
    else
        error('upU-Net not found.')
    end
end

% Apply the network to the test dataset
Y =  predict(net,imdsTest);

for i =1:size(Y,4)
    A = imread(imdsTest.Files{i});
    M = imread(pxdsTest.Files{i});
    montage({double(A),double(M),Y(:,:,1,i)},'DisplayRange',[0,255],'Size',[1,3],...
        'BorderSize',[2,2]);
    colormap parula
    drawnow
    pause

end
