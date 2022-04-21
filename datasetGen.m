%% Dataset generation for upU-Net test
% 
% DATASETGEN applies the procedure described in Benfenati A., upU-Net
% Approaches for Background Emission Remvoal in Fluorescence Microscopy,
% 2022.
%
clearvars
close all
clc

% Select the images size
imgSize = 128;

for i = 1:500
    fprintf('Generating data %3d\n',i);
    [C,B,Cinfo,CLEAR] = particlesSIM(10+randi(5),[76,76,5],[imgSize,imgSize,1],2,...
        'PSF',fspecial('Gaussian',imgSize,2),...
        'sigma',0.03,...
        'poisson',1,...
        'nindex',i);
    % Saving the images in the directories
    imwrite(uint8(B),sprintf('Images/%d/noisy%03d.tiff',imgSize,i));
    imwrite(uint8(CLEAR),sprintf('Labels/%d/noisy%03d.tiff',imgSize,i));
end
