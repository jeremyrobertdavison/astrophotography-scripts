
% ============================================================
% NIGHT SKY STAR COUNTER
% GNU Octave - Astrophotography Image Analysis
%
% Inputs: JPG, JPEG, PNG, TIFF, or BMP image
%
% Outputs:
%   <image_name>_stars_detected.jpg
%   <image_name>_star_coordinates.csv
%
% Requires GNU Octave's image package
% ============================================================

clear;
clc;
close all;

% ------------------------------------------------------------
% 1. SETTINGS
% ------------------------------------------------------------

% 3.0 = More sensitive (may detect more noise)
% 5.0 = Balanced
% 8.0 = Conservative (bright stars only)
sensitivity = 5.0;

% Minimum separation between star peaks (pixels)
min_distance = 5;

% Radius of red circles around detected stars
circle_radius = 9;

% Minimum brightness above local sky background
min_brightness = 0.012;

% ------------------------------------------------------------
% 2. LOAD IMAGE PACKAGE
% ------------------------------------------------------------

try
    pkg load image;
catch
    error(['The Octave image package is required. ', ...
           'Install with: pkg install -forge image']);
end

% ------------------------------------------------------------
% 3. ASK FOR IMAGE FILENAME
% ------------------------------------------------------------

fprintf('\n====================================\n');
fprintf('       NIGHT SKY STAR COUNTER\n');
fprintf('====================================\n\n');

filename = strtrim(input( ...
    'Enter image filename (example: a.jpg): ', 's'));

if isempty(filename)
    error('No filename was entered.');
end

% ------------------------------------------------------------
% 4. LOCATE IMAGE
% ------------------------------------------------------------

% Get folder containing this script
script_path = mfilename('fullpath');

if isempty(script_path)
    script_dir = pwd;
else
    script_dir = fileparts(script_path);
end

% Try the entered filename or full path first
image_path = filename;

if ~exist(image_path, 'file')
    % Try script directory
    image_path = fullfile(script_dir, filename);
end

if ~exist(image_path, 'file')
    error('Image not found: %s', filename);
end

fprintf('\nLoading image: %s\n', image_path);

% Get filename without extension for output files
[image_dir, image_name, image_ext] = fileparts(image_path);

if isempty(image_dir)
    image_dir = pwd;
end

original = imread(image_path);

% ------------------------------------------------------------
% 5. CONVERT IMAGE TO GRAYSCALE
% ------------------------------------------------------------

fprintf('Converting image to grayscale...\n');

if ndims(original) == 3

    gray = 0.2989 * im2double(original(:,:,1)) ...
         + 0.5870 * im2double(original(:,:,2)) ...
         + 0.1140 * im2double(original(:,:,3));

else

    gray = im2double(original);

end

[height, width] = size(gray);

fprintf('Image dimensions: %d x %d pixels\n', ...
        width, height);

% ------------------------------------------------------------
% 6. REDUCE SENSOR NOISE
% ------------------------------------------------------------

fprintf('Reducing image noise...\n');

% Create a small Gaussian blur kernel.
% The kernel suppresses isolated noise while preserving
% many star-like point sources.

[xg, yg] = meshgrid(-2:2, -2:2);

gaussian = exp(-(xg.^2 + yg.^2) / (2 * 0.8^2));

gaussian = gaussian / sum(gaussian(:));

smooth_image = conv2(gray, gaussian, 'same');

% ------------------------------------------------------------
% 7. ESTIMATE BACKGROUND SKY BRIGHTNESS
% ------------------------------------------------------------

fprintf('Estimating sky background...\n');

% Morphological opening removes small bright details
% and approximates the background sky brightness.

background = imopen(smooth_image, ones(31,31));

% Remove estimated sky background
corrected = smooth_image - background;

corrected = max(corrected, 0);

% ------------------------------------------------------------
% 8. ESTIMATE IMAGE NOISE
% ------------------------------------------------------------

fprintf('Calculating detection threshold...\n');

pixels = corrected(:);

median_value = median(pixels);

% Median absolute deviation
mad_value = median(abs(pixels - median_value));

noise_sigma = 1.4826 * mad_value;

% Prevent a zero noise threshold
noise_sigma = max(noise_sigma, 0.0005);

threshold = max(sensitivity * noise_sigma, ...
                min_brightness);

fprintf('Noise estimate: %.6f\n', noise_sigma);

fprintf('Detection threshold: %.6f\n', threshold);

% ------------------------------------------------------------
% 9. IDENTIFY POTENTIAL STARS
% ------------------------------------------------------------

fprintf('Searching for stars...\n');

% Find local maximum brightness in a neighborhood
neighborhood = ones(2 * min_distance + 1);

local_maximum = imdilate(corrected, neighborhood);

% Bright pixels which are local maxima
star_candidates = ...
    (corrected >= local_maximum - 1e-10) ...
    & (corrected > threshold);

% Ignore image boundaries
margin = min_distance + 2;

if height <= 2 * margin || width <= 2 * margin
    error('Image is too small for these detection settings.');
end

star_candidates(1:margin,:) = false;
star_candidates(end-margin+1:end,:) = false;

star_candidates(:,1:margin) = false;
star_candidates(:,end-margin+1:end) = false;

% ------------------------------------------------------------
% 10. GROUP PEAKS INTO STAR DETECTIONS
% ------------------------------------------------------------

fprintf('Validating detected stars...\n');

components = bwconncomp(star_candidates, 8);

% Calculate local 3x3 brightness
aperture_flux = conv2(corrected, ones(3,3), 'same');

star_x = [];
star_y = [];
star_brightness = [];

for k = 1:components.NumObjects

    indices = components.PixelIdxList{k};

    % Locate brightest pixel in this peak
    [~, best] = max(corrected(indices));

    peak_index = indices(best);

    [y, x] = ind2sub(size(gray), peak_index);

    peak_brightness = corrected(y,x);

    % Reject weak detections
    if aperture_flux(y,x) < ...
       max(3 * threshold, 9 * noise_sigma)

        continue;

    end

    % Store detected star
    star_x(end+1,1) = x;
    star_y(end+1,1) = y;

    star_brightness(end+1,1) = peak_brightness;

end

% ------------------------------------------------------------
% 11. SORT DETECTED STARS BY BRIGHTNESS
% ------------------------------------------------------------

[star_brightness, order] = ...
    sort(star_brightness, 'descend');

star_x = star_x(order);
star_y = star_y(order);

star_count = length(star_x);

% ------------------------------------------------------------
% 12. DISPLAY STAR COUNT
% ------------------------------------------------------------

fprintf('\n====================================\n');
fprintf('           ANALYSIS RESULTS\n');
fprintf('====================================\n');

fprintf('Image: %s%s\n', image_name, image_ext);

fprintf('Resolution: %d x %d\n', width, height);

fprintf('Detection sensitivity: %.1f\n', sensitivity);

fprintf('\nTOTAL DETECTED STARS: %d\n', star_count);

fprintf('====================================\n');

% ------------------------------------------------------------
% 13. CREATE ANNOTATED IMAGE
% ------------------------------------------------------------

fprintf('\nMarking stars on image...\n');

% Convert image to RGB 8-bit output
if ndims(original) == 2

    annotated = repmat( ...
        im2uint8(im2double(original)), [1 1 3]);

else

    annotated = im2uint8(im2double(original));

end

% Generate circle coordinates
angles = linspace(0, 2*pi, 72);

% Draw red rings around detected stars
for k = 1:star_count

    cx = star_x(k);
    cy = star_y(k);

    for t = 1:length(angles)

        px = round(cx + ...
             circle_radius * cos(angles(t)));

        py = round(cy + ...
             circle_radius * sin(angles(t)));

        if px >= 1 && px <= width && ...
           py >= 1 && py <= height

            annotated(py,px,1) = 255;
            annotated(py,px,2) = 30;
            annotated(py,px,3) = 30;

        end

    end

end

% ------------------------------------------------------------
% 14. SAVE ANNOTATED IMAGE
% ------------------------------------------------------------

% Save results beside the input image
output_image = fullfile( ...
    image_dir, ...
    [image_name '_stars_detected.jpg']);

imwrite(annotated, output_image, 'Quality', 95);

fprintf('Saved annotated image:\n%s\n', output_image);

% ------------------------------------------------------------
% 15. EXPORT STAR COORDINATES TO CSV
% ------------------------------------------------------------

csv_path = fullfile( ...
    image_dir, ...
    [image_name '_star_coordinates.csv']);

fid = fopen(csv_path, 'w');

if fid == -1
    error('Unable to create CSV file: %s', csv_path);
end

fprintf(fid, 'Star_ID,X_Pixel,Y_Pixel,Brightness\n');

for k = 1:star_count

    fprintf(fid, '%d,%d,%d,%.8f\n', ...
            k, ...
            star_x(k), ...
            star_y(k), ...
            star_brightness(k));

end

fclose(fid);

fprintf('Saved star coordinates:\n%s\n', csv_path);

% ------------------------------------------------------------
% 16. DISPLAY ANNOTATED IMAGE
% ------------------------------------------------------------

fprintf('\nDisplaying image...\n');

figure( ...
    'Name', 'Night Sky Star Detection', ...
    'NumberTitle', 'off');

imshow(annotated);

title(sprintf('%s - Detected Stars: %d', ...
              image_name, star_count));

fprintf('\nAnalysis complete!\n');
fprintf('====================================\n');

