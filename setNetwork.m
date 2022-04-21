function lgraph=setNetwork(imsize,nr_Layers,nr_FirstFilters)
%
% SETNETWORK generates the layer graph for an upU-Net architecture.
%
% graph=setNetwork(imsize,nr_Layers,nr_FirstFilters)
%
% MANDATORY INPUT
% 
% imsize (integer array)    : 2D array with the size of the input images.
%
% nr_Layers (integer)       : number of blocks of layers in the contractive
%                             and expansive parts.
%
% nr_FirstFilters (integer) : number of filters in the very first block.
%                             This number will doubled in the subsequent
%                             blocks.
%
% OUTPUT
%
% lgraph (layer graph)      : layer graph to be trained.
%
%==========================================================================
%
% Version : 1.0 (21-04-2022) Author  : A. Benfenati
% (alessandro.benfenati@unimi.it)
%
%==========================================================================
%
% COPYRIGHT NOTIFICATION
%
% Permission to copy and modify this software and its documentation for
% internal research use is granted, provided that this notice is retained
% thereon and on all copies or modifications. The authors and their
% respective Universities makes no representations as to the suitability
% and operability of this software for any purpose. It is provided "as is"
% without express or implied warranty. Use of this software for commercial
% purposes is expressly prohibited without contacting the authors.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, either visite http://www.gnu.org/licenses/ or
% write to Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
% 02139, USA.
%
%%

% Input layer
convlayers = imageInputLayer(imsize,'Name','input',...
    'Normalization','zscore');


% Contractive path
for i = 1:nr_Layers
    convlayers = [convlayers

    convolution2dLayer(3,nr_FirstFilters,...
    'Name',sprintf('Conv%d',i),...
    'Stride',[2,2],...
    'Padding',[1 1],...
    'PaddingValue','symmetric-exclude-edge');

    batchNormalizationLayer;
    
    reluLayer('Name',sprintf('relu%d',i ))

    ];
    nr_FirstFilters= nr_FirstFilters*2;
end

% Expansive path
for i = 1:nr_Layers
    nr_FirstFilters= nr_FirstFilters/2;
    convlayers = [convlayers

    transposedConv2dLayer(2,nr_FirstFilters,...
    'Name',sprintf('TransConv%d',i),...
    'Stride',2)
    reluLayer('Name',sprintf('relu%d',nr_Layers+i))
    additionLayer(2,'Name',sprintf('add%d',nr_Layers+i));

    ];

end

% Final layers and Regression layer for training
convlayers = [convlayers,
    convolution2dLayer(1,1,'Name','1d1');
    reluLayer('Name','relufin')
    regressionLayer('Name','RegressionOutput')];


% Connecting the contractive path and the expansive path
lgraph = layerGraph(convlayers);
for i=1:nr_Layers
    skipConv = transposedConv2dLayer(2,nr_FirstFilters,...
        'Stride',2,...
        'Name',sprintf('skipConv%d',i));
    nr_FirstFilters= nr_FirstFilters*2;
    lgraph = addLayers(lgraph,skipConv);
end


for i = 1:nr_Layers
    lgraph = connectLayers(lgraph,sprintf('relu%d',i),...
                                  sprintf('skipConv%d',i));
    lgraph = connectLayers(lgraph,sprintf('skipConv%d',i),...
        sprintf('add%d/in2',2*nr_Layers-i+1));
end




