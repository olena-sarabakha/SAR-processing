function [filtered_image] = multilook_filter_sliding(window_dimension,input_image, plot)


%%=========================================================================
% HELP: multilook_filter_sliding function. This function is a multilook filter
% using sliding square window, whose dimensions are defined by the input
% 'window_dimension'. The input image is supposed to be revealed in
% intensity (power). For input variable plot = 0, the function produces an
% image in dB of the output averaged image as well. The sliding window is
% implemented through a 2D convolution.
%%=========================================================================


%% these two lines do the same the two nested cycles do, but muuuuch faster.
window = ones(window_dimension)/(window_dimension^2);
filtered_image = conv2(input_image , window, 'same');


% [rows, columns] = size(input_image); % dimensions of the input image
% 
% % initialization of the counters
% r = 1; % rows counter
% c = 1; % columns counter
% 
% % pre-allocation
% filtered_image = zeros(rows-window_dimension+1, columns-window_dimension+1);
% ratio_image = zeros(rows-window_dimension+1, columns-window_dimension+1);
% 
% % nested cycle in order to have a sliding mechanism
% 
% for i = 1:rows-window_dimension+1 % for cycle on the rows
%     for j = 1:columns-window_dimension+1 % for cycle on the columns 
%         
%         % estimated value of RCS
%         sigma_est = mean(input_image(i:i+ceil(window_dimension/2), j:j+ceil(window_dimension/2)), 'all'); 
%         
%         % construction of the filtered image
%         filtered_image(r, c) = sigma_est; 
%         
%         % ratio image
%         ratio_image(r, c) = input_image(r, c)/filtered_image(r, c); 
%         
%         % counter update
%         c = c+1;
%     end
%     
%     % counter update
%     c = 1;
%     r = r+1;
% end

if plot == 0
    string = sprintf('Moving average - Window %dx%d', window_dimension, window_dimension);
    figure('units','normalized','outerposition',[0 0 1 1]),
    imagesc(10*log10(filtered_image/max(max(filtered_image))));colormap(gray);title(string);
end

end

