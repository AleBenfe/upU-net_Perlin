function [C,B,Cinfo,CLEAR] = particlesSIM(n,ranges,N,r,varargin)
% PARTICLESSIM generates n particles of radius r in the volume of
% dimension ranges(1) X ranges(2) X ranges(3).
%
% [C,B,Cinfo,CLEAR] = particlesSIM(n,ranges,N,r,varargin)
%
% MANDATORY INPUT
% n          : (integer) number of particles
% ranges     : (double array) 3D array containg the dimension of the volume
% N          : (double arrayr) 3D array containg the number of voxel in
%                              each dimension
% r          : (double) radius of the particles
%
% OPTIONAL INPUT
%It needs improvement
% SIGMA     : (double)       standard deviation of the Gaussian noise
%                            affecting the data.
%                            Default: 0;
% POISSON   : (double)       flag indicating if Poisson noise affects data.
%                            Default: 0
% PSF       : (double)       PSF for the blurr.
%                            Default: Gaussian.
% CENTER    : (double array) initial position of the beads.
%                            Default: random.
%
% OUTPUT
%
% C         : (double array) Center of the particles
% Cinfo     : (structure)    structure containing some infos for
%                            performance evaluation
% B         : (double array) 3D array containing the images stacked
%                            vertically
% CLEAR     : (double array) 3D array containing the clean images stacked
%                            vertically
%
%==========================================================================
%
% Version : 1.0 (21-04-2022)
% Author  : A. Benfenati (alessandro.benfenati@unimi.it)
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


% Unpack data
rangeX  = ranges(1);
rangeY  = ranges(2);
rangeZ  = ranges(3);
Nx      = N(1);
Ny      = N(2);
Nz      = N(3);
dx      = rangeX/Nx;
dy      = rangeY/Ny;
dz      = rangeZ/Nz;


%% Defaults
sigma       = 0;
poisson     = 0;
C_init      = rand(n,3)*diag(ranges);
filt        = fspecial('gaussian',5,1) ;

if (nargin-length(varargin)) ~= 4
    error('Wrong number of required parameters');
end

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
    fprintf('!\n');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SIGMA'
                sigma = varargin{i+1};
            case 'PSF'
                filt = varargin{i+1};
            case 'POISSON'
                poisson = varargin{i+1};
            case 'CENTER'
                C_init = varargin{i+1};
        end
    end
end


if Nx~=Ny
    error('Not handled yet... Sorry')
end

% Data creation
B      = zeros(N(1),N(2),N(3));

% Remove intersecting particles
D      = (pdist(C_init));
while min(D(:))<2*r 
    mnm = min(D(:));
    D = squareform(D);
    [i,~] = find(D==mnm);
    C_init(i(1),:) = rand(1,3)*diag(ranges);
    D = (pdist(C_init));
end



C      = C_init;

% Benchmark for control
Cinfo = cell(n,1);
for k=1:n
    Cinfo{k} = struct('Pos',[],'SpanFrames',[],'Radii',[]);
end


% Coordinates for circle discretization
[X,Y]   = meshgrid(1:N(1),1:N(2));
X       = dx*X(:);
Y       = dy*Y(:);
Z       = kron((1:N(3))',dz*ones(N(1)*N(2),1));
data    = [kron(ones(N(3),1),X),kron(ones(N(3),1),Y), Z ];

% Geenrate Perlin noise
A = 0.9*uint8(255*perlinNoise3D([Nx, Ny]));


CLEAR   = zeros(prod(N),1);
V       = [kron(ones(N(3),1),A(:))];

for k = 1:n
    dist    = sum((data-kron(C(k,:),ones(size(data,1),1))).^2,2)<r^2;
    % 0.9*max(A(:));
    V(dist)     = 229.5;
    CLEAR(dist) = 0.9*max(A(:));
end

% Reshape the data in the correct dimension
V = reshape(double(V),N);

CLEAR   = reshape(CLEAR,N);
V       = imfilter(V,filt,'symmetric');
CLEAR   = imfilter(CLEAR,filt,'symmetric');

for z = 1:Nz
    for k = 1:n
        % Check the intersection of the k--th sphere with the z--th plane
        l = abs(C(k,end)-z*dz);
        if (l<r && abs(l-r)>1e-8)
            % Save the useful information about this center
            Cinfo{k}.Pos          = C(k,1:2);
            Cinfo{k}.SpanFrames   = [Cinfo{k}.SpanFrames; z];
            Cinfo{k}.Radii        = [Cinfo{k}.Radii; sqrt(r^2-l^2)];
        end
    end
    % Gaussian Noise
    if sigma>0
        noise    = randn(Nx,Nx);
        noise    = noise/norm(noise,'fro')*(norm(V(:,:,z),'fro'));
        V(:,:,z) = V(:,:,z) + sigma*noise;
        V(:,:,z) = V(:,:,z) + abs(min(min((V(:,:,z)))));
        CLEAR = CLEAR + sigma*noise;
    end
end
% Poisson noise
if poisson
    V = 1e12*imnoise(1e-12*V,'Poisson'); 
    %CLEAR = 1e12*imnoise(1e-12*CLEAR,'Poisson'); 
end



% srhitnking the value in the range [0,255]
CLEAR = (CLEAR-min(CLEAR(:)))/(max(CLEAR(:))-min(CLEAR(:)));
CLEAR = uint8(CLEAR*255);

V = (V-min(V(:)))/(max(V(:))-min(V(:)));
V = V*255;

B(:,:,:) = uint8(V);