
%%=========================================================================
%%HELP: function to resample an image. The new sample numbers are chosen in
%%such a way to be a power of 2, in order to make the FFT faster. Inputs:
%%image (in double time domain), multiplier (multipliying factor to
%%obtain the new number of samples; a suggested number is 4, or any power
%%of 2). Outputs: resampled_image in double time domain, rows_resampled
%%(new number of range samples), columns_resampled (new number of azimuth
%%samples). If multiplier_rows or multiplier_columns is 0, the matrix won't
%%be resampled in that dimension.

%%tau = fast time; ta = slow time; fr = range frequency; fa = azimuth
%%frequency.

%%=========================================================================

function [resampled_image, rows_resampled, columns_resampled]=resampling(image, multiplier_rows, multiplier_columns)

[rows, columns] = size(image);

% if the number of inputs is 2, the multipliers of rows and columns are
% assumed to be the same
if nargin == 2
   multiplier_columns = multiplier_rows;
end

% definition of the new number of rows
rows_resampled = multiplier_rows*(2^nextpow2(rows)); % new number of samples in range
if multiplier_rows == 0 
    rows_resampled = rows;
end
% definition of the new number of columns
columns_resampled = multiplier_columns*(2^nextpow2(columns)); % new number of samples in azimuth
if multiplier_columns == 0
    columns_resampled = columns;
end

% range padding
image_fr_ta = fftshift(fft(ifftshift(image, 1), [], 1), 1); % FFT on rows (range) -->(fr, ta)
image_fr_ta = [zeros((rows_resampled-rows)/2, columns); image_fr_ta; zeros((rows_resampled-rows)/2, columns)]; % zero padding
resampled_image_aux = ifftshift(ifft(ifftshift(image_fr_ta, 1), [], 1), 1); % IFFT on rows (range) --> (tau, ta)

% azimuth padding
image_tau_fa = fftshift(fft(ifftshift(resampled_image_aux, 2), [], 2), 2); % FFT on columns (azimuth) --> (tau, fa)
image_tau_fa = [zeros(rows_resampled, (columns_resampled-columns)/2) image_tau_fa zeros(rows_resampled, (columns_resampled-columns)/2)]; % zero padding
resampled_image = ifftshift(ifft(ifftshift(image_tau_fa, 2), [], 2), 2); % IFFT on columns (azimuth) --> (tau, ta)

end