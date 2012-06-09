% Hacked together code for implementing the colour magnifying component
% of Eulerian Video Magnification:
%   http://people.csail.mit.edu/mrub/vidmag/
%
% Only enhances the red RGB component at the moment. It can be
% extended to handle all three; memory limitations on my laptop
% prevents me from doing them all in one step.
%
% Requires Minh N. Do's laplacian pyramid toolbox:
%   http://www.ifp.illinois.edu/~minhdo/software/lptoolbox.tar.gz
%
%   Copyright 2012
%       James Brotchie <brotchie@gmail.com>
%       https://github.com/brotchie/
%

HOME = getenv('HOME');
LAPLACIAN_PYRAMID_TOOLBOX_PATH = [HOME filesep 'Downloads' filesep 'image'];
addpath(LAPLACIAN_PYRAMID_TOOLBOX_PATH);

% Video to process.
VIDEO_PATH = [HOME filesep 'Downloads' filesep 'IMG_0940.MOV'];

% Depth of the laplacian pyramid.
NDECOMPOSITIONS = 4;

% Scaling factor for bandpass filtered signal.
SCALING_FACTOR_ALPHA = 5;

% Type of Laplacian decomposition.
DECOMPOSITION_TYPE = '9-7';

% Filter setup.
FILTER_TYPE = 'butter';
FILTER_ORDER = 32;
FILTER_LP_CUTOFF_HZ = 0.5;
FILTER_HP_CUTOFF_HZ = 2;

video = VideoReader(VIDEO_PATH);

% Read video frame, we permute to flip from landscape
% to portrait.
frames = permute(video.read(), [2 1 3 4]);

VIDEO_HEIGHT = size(frames, 1);
VIDEO_WIDTH = size(frames, 2);

assert(VIDEO_HEIGHT > VIDEO_WIDTH, ...
       'We assume that the video width is greater than the height when padding.');

% Pad the video so it's square; laplacian pyramid toolbox requires
% square inputs.
PADDING = (VIDEO_HEIGHT - VIDEO_WIDTH) / 2;
frames = [zeros(VIDEO_HEIGHT, PADDING, 3, size(frames,4), 'uint8') ...
          frames ...
          zeros(VIDEO_HEIGHT, PADDING , 3, size(frames,4), 'uint8')];
nFrames = size(frames, 4);

% Pre-allocate memory to store red channel decompositions.
decomposed = cell(NDECOMPOSITIONS+1, 1);
for ii = 1:NDECOMPOSITIONS+1
    N = size(frames, 1) / 2^(NDECOMPOSITIONS - ii + 1);
    decomposed{ii} = zeros(N, N, nFrames);
end

% Perform the decompositions and store them in the decomposed cell
% array.
for ii = 1:nFrames
    fprintf('Decomposing frame %d of %d.\n', ii, nFrames);
    y = lpd(double(frames(:, :, 1, ii)), DECOMPOSITION_TYPE, NDECOMPOSITIONS);
    for jj = 1:NDECOMPOSITIONS
        decomposed{jj}(:,:,ii) = y{jj};
    end
end

% Filter each decomposed pixel in the time domain.
f = fdesign.bandpass('n,fc1,fc2', FILTER_ORDER, ...
                     FILTER_LP_CUTOFF_HZ, FILTER_HP_CUTOFF_HZ, ...
                     video.FrameRate);
h = design(f, FILTER_TYPE);

improved = cell(NDECOMPOSITIONS+1, 1);
for ii = 1:NDECOMPOSITIONS+1
    fprintf('Band pass filtering laplacian pyramid level %d of %d.\n', ii, NDECOMPOSITIONS+1);
    improved{ii} = decomposed{ii} + SCALING_FACTOR_ALPHA*filter(h, decomposed{ii}, 3);
end

% Clear the decomposed frames to free up a bit of memory.
clear decomposed;

% Recompose the filtered signal.
red_signal = zeros(480, 480, nFrames, 'uint8');
M = repmat(struct('cdata', zeros(size(red_signal, 1), size(red_signal, 2), 3, 'uint8'), ...
                  'colormap', []), [nFrames 1]);

for ii = 1:nFrames
    fprintf('Recomposing frame %d of %d.\n', ii, nFrames);
    y = cellfun(@(x)x(:,:,ii), improved, 'UniformOutput', 0);
    signal(:,:,ii) = uint8(lpr(y, DECOMPOSITION_TYPE));
    M(ii).cdata = cat(3, signal(:,:,ii),frames(:, :, 2:3, ii));
end

figure;
movie(M, 1, video.FrameRate);
