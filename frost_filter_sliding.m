function [filtered_image, ratio_image] = frost_filter_sliding(window_dimension, input_image, ENL, plot)

%%=========================================================================
%HELP: function that applies the Forst filter with a sliding window. Input:
%window_dimension, dimension of the square window; input_image, image to
%despeckle, revealed in intensity; ENL, equivalent number of looks; plot,
%if 0 an image of the output and ratio matrices generated. Output:
%filtered_image, image filtered by the Frost filter; ratio_image, ratio
%between input and output image.
%%=========================================================================

[rows, columns] = size(input_image); % dimensions of the input image

% initialization of the counters
r = 1; % rows counter
c = 1; % columns counter

% pre-allocation
filtered_image = zeros(rows-window_dimension+1, columns-window_dimension+1);
ratio_image = zeros(rows-window_dimension+1, columns-window_dimension+1);

% nested cycle in order to have a sliding mechanism

for i = 1:rows-window_dimension+1 % for cycle on the rows
    
    for j = 1:columns-window_dimension+1 % for cycle on the columns
        I_local = input_image(i:i+window_dimension-1, j:j+window_dimension-1); % selected window of the image
        local_variance = var(I_local, 1, 'all'); % variance of the selected window
        local_mean = mean(I_local, 'all'); % mean value of the selected window
        local_contrast = local_variance/local_mean^2; % local contrast
        
        % calculation of k parameter
        k = (local_contrast-1/ENL)/(local_contrast*(1+1/ENL)); 
        % Controllo sul segno di k
        if k < 0
            k = 0;
        end
        
        % estimated value of RCS
        sigma_est = local_mean+k*(input_image(i+ceil(window_dimension/2), j+ceil(window_dimension/2))-local_mean); 
        
        % construction of the filtered image
        filtered_image(r, c) = sigma_est; 
        
        % ratio image
        ratio_image(r, c) = input_image(r, c)/filtered_image(r, c); 
        
        % counter update
        c = c+1;
    end
    
    % counter update
    c = 1;
    r = r+1;
end

% generation of plot if required
if plot == 0
    string = sprintf('Sliding Frost filter - Window %dx%d', window_dimension, window_dimension);
    figure('units','normalized','outerposition',[0 0 1 1]),
    imagesc(10*log10(filtered_image/max(max(filtered_image))));colormap(gray);title(string);
    
    string = sprintf('Ratio image from sliding Frost filter - Window %dx%d', window_dimension, window_dimension);
    figure('units','normalized','outerposition',[0 0 1 1]),
    imagesc(10*log10(ratio_image/max(max(ratio_image))));colormap(gray);title(string);
end


end
        