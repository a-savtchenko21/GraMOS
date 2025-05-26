
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Copyright Disclaimer                                    %
%               Written by Michael Reiss                                %
%               m1reiss@ucsd.edu                                        %
%               University of California San Diego                      %
%               Department of Bioengineering                            %
%               MIT license                                             %
%               Last date modified: March 6, 2025                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;

%% Read in data
%
data = single(tiffreadVolume('/Users/michaelreiss/Documents/UCSD/FaceMaps/Termination Project/InVivo/CalciumSignaling/2025-01-20_WTAD st_1.tif')); % CHANGE THIS PATH TO THE INPUT VIDEO TO BE ANALYZED
center      = [533, 556];   % WT
start_stim = 64; % WT, Frame in which the stimulus starts 
min_clim_factor = 4; % For plotting
max_clim_factor = 0.07; % For plotting
Fs = 25;
%}

%{
data = single(tiffreadVolume('/Users/michaelreiss/Documents/UCSD/FaceMaps/Termination Project/InVivo/CalciumSignaling/2025-01-20_M233L st.tif'));
center      = [523,533];  
start_stim  = 165;          
min_clim_factor = 2.5;
max_clim_factor = 0.15;
Fs = 25;
%}

%{
data = single(tiffreadVolume('/Users/michaelreiss/Documents/UCSD/FaceMaps/Termination Project/InVivo/CalciumSignaling/new videos/2025-01-20_WTAD.tif'));
center      = [534,558]; 
start_stim  = 98;          
min_clim_factor = 100; % min value is 0, so large factor needed
max_clim_factor = 0.3;
Fs = 33.5;
%}

%{
data = single(tiffreadVolume('/Users/michaelreiss/Documents/UCSD/FaceMaps/Termination Project/InVivo/CalciumSignaling/new videos/2025-01-20_M233L st_19.tif'));
center      = [531,555];  
start_stim  = 109;          
min_clim_factor = 4;
max_clim_factor = 0.15;
Fs = 25;
%}

%{
data = single(tiffreadVolume('/Users/michaelreiss/Documents/UCSD/FaceMaps/Termination Project/InVivo/CalciumSignaling/Paper/dataset3/2025-01-20_M233L st_22.tif'));
center      = [526,555];  
start_stim  = 48;          
min_clim_factor = 3.5;
max_clim_factor = 0.05;
Fs = 25;
%}

%{
data = single(tiffreadVolume('/Users/michaelreiss/Documents/UCSD/FaceMaps/Termination Project/InVivo/CalciumSignaling/Paper/dataset3/2025-01-20_WTAD st_7.tif'));
center      = [530,555];  
start_stim  = 77;          
min_clim_factor = 100;
max_clim_factor = 0.3;
Fs = 25;
%}



% Get data size
n_rows = size(data,1);
n_colums = size(data,2);
n_samples = size(data,3);



data_img = data; % Data as image sequences

%% Filter data
[U,S,V]=svd(reshape(data,[n_rows*n_colums,n_samples]),0);
        sig_modes = 20; % Choose eigenmodes depending on eigenspectrum. Analysis is not very sensitive to that
        Un=single(U(:,1:sig_modes));
        Sn=single(S(1:sig_modes,1:sig_modes));
        Vn=single(V(:,1:sig_modes));
data = Un*Sn*Vn';
clear U V S Un  Sn Vn


%% Choose only regions that show some activity -> reduces computational time
max_values = max(data');
min_values = min(data');

abs_diff = max_values-min_values;
no_plot = abs_diff<20; % Can be changed depending on data
toplot = abs_diff>=20; % Can be changed depending on data

data_ROI = data(toplot,:); % Select only data of interest


%% Get biggest activation out
spikevisuals = zeros(size(data_ROI,1), n_samples);
for i = 1:size(data_ROI,1)
    % Calculate the angle of the hilbert transform of the smoothened data
    phase_calc = phase(hilbert(smooth(data_ROI(i,:),40)));
    % Spiketimings are minima in the angle that are here at least 20 frames
    % away and have a prominence of at least 0.05
    spikes = find(islocalmin(phase_calc,'MinSeparation', 20, 'MinProminence', 0.05));
    spikes(spikes>n_samples-30) = []; % Exclude very late spikes that are artifacts of the hilbert transform and angular calucalation
    spikes(spikes<30) = [];
    % If only 1 spike (the main spike from the stimulus) is wanted,
    % un-comment
    %{
    [~, phase_minidx] = min(phase_calc);
    spikes = intersect(spikes, phase_minidx);
    %}

    spikevisuals(i,spikes) = 1;
end


%% Plot number of spikes per pixel
%
spikes_per_pixel_ROI = sum(spikevisuals,2);
spikesmap = zeros(1,n_rows*n_colums);
spikesmap(toplot) = spikes_per_pixel_ROI;
spikesmap = reshape(spikesmap, [n_rows, n_colums]);
spikesmap(spikesmap>4) = 4; % Limit to maximal 4 spikes

min_clim = min_clim_factor*(min(data_img(:))+1); 
max_clim = max_clim_factor*max(data_img(:)); %
plot_background = squeeze(data_img(:,:,50)); % choose frame 50 for background
plot_foreground = spikesmap;

figure(1);
ax1 = axes;
imagesc(ax1, plot_background);
axis equal
axis tight
colormap(ax1, 'abyss');     % background colormap
clim(ax1, [min_clim, max_clim]);  
%title(ax1, 'Background');
hold(ax1, 'on');

%% Axes #2 for overlay
ax2 = axes;
imagesc(ax2, plot_foreground);
axis equal
axis tight
colormap(ax2, 'jet');      % overlay colormap
clim(ax2, [0, 4]);        % color limit for data2

% Make the background of ax2 transparent
set(ax2, 'Color', 'none');
% Turn off the axis lines for ax2 (optional, so you see only one set of ticks)
set(ax2, 'XColor', 'none', 'YColor', 'none');

% Make the zero-values transparent in the overlay:
%   h2 is the handle to the imagesc in ax2
h2 = ax2.Children;  % the imagesc is typically ax2.Children(1)
set(h2, 'AlphaData', plot_foreground ~= 0);

% Link the two axes so they stay in sync (zoom, pan, etc.)
linkaxes([ax1, ax2], 'xy');

% Define scale bar parameters
scale_bar_um = 100;             % scale bar length in µm
scale_bar_pixels = scale_bar_um / 1.2;  % convert µm to pixels (~83 pixels)
% Define the position of the scale bar (adjust margins as needed)
margin = 20;  % margin in pixels from the image border
x_start = margin;  
y_start = size(plot_background,1) - margin;  % bottom of image (in pixel coordinates)
% Plot the scale bar on the background axes (ax1)
hold(ax1, 'on');
line(ax1, [x_start, x_start + scale_bar_pixels], [y_start, y_start], 'Color', 'w', 'LineWidth', 3);
% Add a text label above the scale bar
text(ax1, x_start, y_start - 15, sprintf('%d \\mum', scale_bar_um), ...
    'Color', 'w', 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
%}

%{
%% Create a tail for each spike for the video
for i = 1:size(data_ROI,1)

    spiketimings = find(spikevisuals(i,:)==1);
    
    %% 3) For each spike, fill subsequent values with a decreasing 0.99, 0.98, ...
    for k = 1:length(spiketimings)
        st = spiketimings(k);  % the time of the spike
        nextIndex = st + 1;
        value = 0.99;          % start the decay
        
        % continue until the end of the time axis or until value reaches 0
        while nextIndex <= n_samples && value > 0
            % Option B (alternative): if you want to keep the maximum
            % in case of overlapping decays, uncomment below and comment above
            % spikevisuals(i, nextIndex) = max(spikevisuals(i, nextIndex), value);
            spikevisuals(i, nextIndex) = max(spikevisuals(i, nextIndex), value);
            % decrement the value and move to the next time index
            value = value - 0.01;
            nextIndex = nextIndex + 1;
        end
    end
end
%}

plot_video = zeros(n_rows*n_colums,n_samples);
plot_video(toplot,:) = spikevisuals;

plot_video = reshape(plot_video, [n_rows,n_colums,n_samples]);


%% Plot video IF DESIRED 
%{
for i = 1:n_samples
    figure(1)
    imagesc(squeeze(plot_video(:,:,i)));
    clim([0,1])
    colormap('jet')
    title(num2str(i))
    %pause(0.1)
    graphics(i) = getframe(gcf); % Use this to save the frames 
    
end
   
v = VideoWriter('check.mp4', 'MPEG-4'); % CHANGE NAME AS DESIRED
v.FrameRate = 20;
%v.Quality = 100;
open(v)

for i = 1:n_samples
    writeVideo(v, graphics(i));
end

close(v);
%}



min_clim = min_clim_factor*(min(data_img(:))+1); 
max_clim = max_clim_factor*max(data_img(:)); 

%{

%% Axes #1 for background

for i = 1:n_samples
    plot_background = squeeze(data_img(:,:,i));
    plot_foreground = squeeze(plot_video(:,:,i));

    figure(1)
    ax1 = axes;
    imagesc(ax1, plot_background);
    axis equal
    axis tight
    colormap(ax1, 'abyss');     % background colormap
    clim(ax1, [min_clim, max_clim]);  
    %title(ax1, 'Background');
    hold(ax1, 'on');
    
    %% Axes #2 for overlay
    ax2 = axes;
    imagesc(ax2, plot_foreground);
    axis equal
    axis tight
    colormap(ax2, 'jet');      % overlay colormap
    clim(ax2, [0 1]);        % color limit for data2
    
    % Make the background of ax2 transparent
    set(ax2, 'Color', 'none');
    % Turn off the axis lines for ax2 (optional, so you see only one set of ticks)
    set(ax2, 'XColor', 'none', 'YColor', 'none');
    
    % Make the zero-values transparent in the overlay:
    %   h2 is the handle to the imagesc in ax2
    h2 = ax2.Children;  % the imagesc is typically ax2.Children(1)
    set(h2, 'AlphaData', plot_foreground ~= 0);
    
    % Link the two axes so they stay in sync (zoom, pan, etc.)
    linkaxes([ax1, ax2], 'xy');

    graphics(i) = getframe(gcf);
end


%}



%% Help calucate the number of pixels/neurons that fire how many times
%{
a = zeros(1,size(data_ROI,1));
for i = 1:size(data_ROI,1)

    spiketimings = find(spikevisuals(i,:)==1);
a(i) = length(spiketimings);
end

a(a==0) = [];

b = a;
b(b<2) = [];


histogram(a)
xlabel('Number Activated pixels')
ylabel('Count')
act_percentage = length(b)/length(a);
%}



%% Calculate aboslute spiketiming of the entire video
%
activation_times = zeros(1, size(data_ROI,1));
for i = 1:size(data_ROI,1)
    phase_calc = phase(hilbert(smooth(data_ROI(i,:),40)));
    spikes = find(islocalmin(phase_calc,'MinSeparation', 20, 'MinProminence', 0.05));
    spikes(spikes>n_samples-30) = []; % Exclude very late spikes that are artifacts of the hilbert transform and angular calucalation
    spikes(spikes<30) = [];

    % If only 1 spike is wanted
    [~, phase_minidx] = min(phase_calc);
    spikes = intersect(spikes, phase_minidx);
    %}
    if ~isempty(spikes)
        activation_times(i) = spikes;
    end
end

activation_map = zeros(1,n_rows*n_colums);
activation_map(toplot) = activation_times;

activation_map = reshape(activation_map, [n_rows,n_colums]);




%% Plot
min_clim = min_clim_factor*(min(data_img(:))+1); 
max_clim = max_clim_factor*max(data_img(:)); 
plot_background = squeeze(data_img(:,:,50)); % choose frame 50 as background
plot_foreground = activation_map./Fs;

figure(2)
ax1 = axes;
imagesc(ax1, plot_background);
axis equal
axis tight
colormap(ax1, 'abyss');     % background colormap
clim(ax1, [min_clim, max_clim]);  
%title(ax1, 'Background');
hold(ax1, 'on');
%% Axes #2 for overlay
ax2 = axes;
imagesc(ax2, plot_foreground);
axis equal
axis tight
colormap(ax2, 'jet');      % overlay colormap
clim(ax2, [min(plot_foreground(:)), max(plot_foreground(:))]);        % color limit for data2
% Make the background of ax2 transparent
set(ax2, 'Color', 'none');
% Turn off the axis lines for ax2 (optional, so you see only one set of ticks)
set(ax2, 'XColor', 'none', 'YColor', 'none');
% Make the zero-values transparent in the overlay:
%   h2 is the handle to the imagesc in ax2
h2 = ax2.Children;  % the imagesc is typically ax2.Children(1)
set(h2, 'AlphaData', plot_foreground ~= 0);
% Link the two axes so they stay in sync (zoom, pan, etc.)
linkaxes([ax1, ax2], 'xy');
% Save the original position of the overlay axes
ax2_pos = get(ax2, 'Position');
% Create the colorbar for the foreground image (associated with ax2)
cb = colorbar(ax2);
% Adjust the colorbar's position so that it sits to the right of ax2 without shifting it.
% Here, we place the colorbar just to the right of ax2 by modifying its Position property.
cb.Position = [ax2_pos(1) + ax2_pos(3) - 0.08, ax2_pos(2), 0.02, ax2_pos(4)];
% Define scale bar parameters
scale_bar_um = 100;             % scale bar length in µm
scale_bar_pixels = scale_bar_um / 1.2;  % convert µm to pixels (~83 pixels)
% Define the position of the scale bar (adjust margins as needed)
margin = 20;  % margin in pixels from the image border
x_start = margin;  
y_start = size(plot_background,1) - margin;  % bottom of image (in pixel coordinates)
% Plot the scale bar on the background axes (ax1)
hold(ax1, 'on');
line(ax1, [x_start, x_start + scale_bar_pixels], [y_start, y_start], 'Color', 'w', 'LineWidth', 3);
% Add a text label above the scale bar
text(ax1, x_start, y_start - 15, sprintf('%d \\mum', scale_bar_um), ...
    'Color', 'w', 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');


%}


%% Calculate gradient of activation map
%{
activation_map = NaN(1,n_rows*n_colums);
activation_map(toplot) = activation_times;

activation_map = reshape(activation_map, [n_rows,n_colums]);
activation_map(activation_map<=start_stim) = NaN;

% Step 1: Compute partial derivatives (gradient) of T wrt rows, cols
[Ty, Tx] = gradient(activation_map); 

% Exclude gradients that definitely do not show connectivity
exclude = find(abs(Tx)>10 | abs(Ty)>10);
Tx(exclude) = NaN;
Ty(exclude) = NaN;


% Note: 'gradient' in MATLAB with [Fx, Fy] = gradient(F) 
% => Fx is ∂F/∂x (difference in columns), 
%    Fy is ∂F/∂y (difference in rows).
% Because row index ~ y, column index ~ x, 
% we have to be careful with the ordering returned by gradient().
%
% Here, I've used [Ty, Tx] just to keep consistent naming:
%   Ty = ∂(T)/∂(row)
%   Tx = ∂(T)/∂(col)

% Step 2: Compute speed = 1 / sqrt((dT/dx)^2 + (dT/dy)^2)
% We'll add a small epsilon to avoid division by zero
epsilon = 1e-12;
speed_map = 1 ./ sqrt(Tx.^2 + Ty.^2 + epsilon);
%limit max speed
speed_map(speed_map>2) = 2;

% Step 3: Where original was NaN, keep it NaN 
% (This may already be the case if the gradient is NaN, 
%  but let's be explicit.)
nanMask = isnan(activation_map);
speed_map(nanMask) = NaN;

% Step 4: Visualize
figure;
imagesc(speed_map); 
axis image;  % so the aspect ratio is 1:1
colormap(jet); 
clim([0,1])
colorbar;
axis tight
axis equal
%}





%% Calculate propagatino speed by selcting activation times 100+-5 pixels away from the center stimulus
%{
activation_map = activation_map-start_stim;
radius = 100;
tolerance = 0.05;     % ±5%
rMin = radius * (1 - tolerance);  
rMax = radius * (1 + tolerance);  

% Distance map
[nRows, nCols] = size(activation_map);
[colIdx, rowIdx] = meshgrid(1:nCols, 1:nRows);
distMap = sqrt( (rowIdx - center(1)).^2 + (colIdx - center(2)).^2 );

% Mask for ring [rMin, rMax]
ringMask = distMap >= rMin & distMap <= rMax;

% Extract values in ring, remove NaNs
vals_ring = activation_map(ringMask);
vals_ring = vals_ring(~isnan(vals_ring));

% Plot histogram
figure;
histogram(vals_ring, 50);
title(sprintf('Histogram of activation timings after stimulus: radius = %d ± %.1f%%', radius, tolerance*100));
xlabel('Activation time (frames)');
ylabel('Count of pixels');



%{
%% Gaussian fitting here
data = vals_ring;
% data is a vector of raw observations
numGauss = 4; % number of gaussians to fit
gmModel = fitgmdist(data, numGauss, 'Options', statset('MaxIter',1000));

figure; 
histogram(data, 'Normalization','pdf');  % plot histogram as pdf
hold on;
xgrid = linspace(min(data), max(data), 200);
pdfVals = pdf(gmModel, xgrid(:)); 
plot(xgrid, pdfVals, 'r-', 'LineWidth',2);
title('Data + 2-Gaussian Mixture PDF');
legend('Data histogram','GMM PDF');

means = gmModel.mu(:); stds = sqrt(squeeze(gmModel.Sigma));
%}
%}