function v = perlinNoise3D(N,varargin)
%
%PERLINNOISE3D creates a 3D array containing Perlin noise
%
%  v= perlinNoise3D(N,varargin)
%
% MANDATORY INPUT
%
% N (integer array) : a array with 3 elements with the size of the desired
%                     output
%
% OPTIONAL INPUT
%
% f (double)        : the initial frequency of the noise. It is doubled on
%                     each octave agfter the first. Default: 1.0.
%
% a (double)        : the initial amplitude of the noise. Default. 1.0
%
% p (double)        : the persistence, i.e. the factor to which the
%                     amplitute is multiplited for for deeper octaves.
%                     Default: 0.5.
%
% x_scale (double)  : scaling in the x-axis. Default: 1.
%
% y_scale (double)  : scaling in the y-axis. Default: 1.
%
% z_scale (double)  : scaling in the z-axis. Default: 1.
%
%
% OUTPUT
% 
% v (double array)  : array of the chosen dimension with Perlin noise with
%                     values in [0,1].
%
%==========================================================================
%
% Version : 1.1 (02-02-2023) Author  : A. Benfenati
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

if max(size(N))==2
    N(3) = 1;
end
if prod(size(N))>3
    error('Dimension mismatch: 4D noise (or larger) not implemented yet.')
end

freqInit    = 1.0;
amplInit    = 1.0;
persistence = 0.5;
x_scale     = 1;
y_scale     = 1;
z_scale     = 1;

if (nargin-length(varargin)) ~= 1
    error('Wrong number of required parameters');
end

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
    fprintf('!\n');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'F'
                freqInit = varargin{i+1};
            case 'A'
                amplInit = varargin{i+1};
            case 'P'
                persistence = varargin{i+1};
            case 'X_SCALE'
                x_scale = varargin{i+1};
            case 'Y_SCALE'
                y_scale = varargin{i+1};
            case 'Z_SCALE'
                z_scale = varargin{i+1};
        end
    end
end


pos = randperm(256);
pos = [pos pos];

v = zeros(N);

for z = 1:N(3)
    for x = 1:N(1)
        for y = 1:N(2)

            temp      = 0;
            frequency = freqInit;
            amplitude = amplInit;

            xf = x_scale*(x-1)/N(1)+1;
            yf = y_scale*(y-1)/N(2)+1;
            zf = z_scale*(z-1)/N(3)+1;
            for k = 1:8
                temp = temp + amplitude*truePerlin3d(frequency*xf,frequency*yf,frequency*zf,pos);
                frequency = frequency*2;
                amplitude = amplitude*persistence;
            end

            v(x,y,z)  = temp;
        end
    end
end

% Data rescaled in [0,1];
v = v-min(v(:));
v = v/max(v(:));
end

function n = truePerlin3d(x,y,z,pos)

X = mod(floor(x),255)+1;
Y = mod(floor(y),255)+1;
Z = mod(floor(z),255)+1;

xi = x - floor(x);
yi = y - floor(y);
zi = z - floor(z);

p000 = pos(pos(pos(X)+Y)+Z);
p100 = pos(pos(pos(X+1)+Y)+Z);
p110 = pos(pos(pos(X+1)+Y+1)+Z);
p010 = pos(pos(pos(X)+Y+1)+Z);
p001 = pos(pos(pos(X)+Y)+Z+1);
p101 = pos(pos(pos(X+1)+Y)+Z+1);
p111 = pos(pos(pos(X+1)+Y+1)+Z+1);
p011 = pos(pos(pos(X)+Y+1)+Z+1);

v000  = dot(getVector(p000),[xi,yi,zi]);
v100  = dot(getVector(p100),[xi-1,yi,zi]);
v110  = dot(getVector(p110),[xi-1,yi-1,zi]);
v010  = dot(getVector(p010),[xi,yi-1,zi]);
v001  = dot(getVector(p001),[xi,yi,zi-1]);
v101  = dot(getVector(p101),[xi-1,yi,zi-1]);
v111  = dot(getVector(p111),[xi-1,yi-1,zi-1]);
v011  = dot(getVector(p011),[xi,yi-1,zi-1]);

nb = lerp(lerp(v000,v100,xi),...
    lerp(v010,v110,xi),...
    yi);
nu = lerp(lerp(v001,v101,xi),...
    lerp(v011,v111,xi),...
    yi);
n = lerp(nb,nu,zi);

end

function g = getVector(x)

h = mod(x,8);

switch h
    case 0
        g = [1 1 1];
    case 1
        g = [-1 1 1];
    case 2
        g = [1 -1 1];
    case 3
        g = [1 1 -1];
    case 4
        g = [-1 -1 1];
    case 5
        g = [-1 1 -1];
    case 6
        g = [1 -1 -1];
    case 7
        g = [-1 -1 -1];
end

end

function y = lerp(a,b,x)
y = (1-fade(x))*a + b*fade(x);
end

function y = fade(x)
y = x * x * x * (x * (x * 6 - 15) + 10);
end

