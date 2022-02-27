
%%=========================================================================
%%HELP: function that compresses the input image in the azimuth domain with
%%a matched filter.
%%Input: raw_image, image that has to be compressed; beta, focalization
%%parameter, calculated as beta=(2*pi*va^2)/(lambda*R0), where va is the
%%platform's velocity, lambda is the wavelength of the signal, R0 is the
%%closest-approach distance; fa, is the frequency vector in azimuth.
%%Output: filtered_image, input image compressed in azimuth.

%%tau is the fast time (range), ta is the slow time (azimuth). The input
%%image is supposed to be in the double time domain (tau, ta).
%%=========================================================================

function [filtered_image]=matched_filter_azimuth(raw_image, beta, fa)

[rows, ~] = size(raw_image); % raw image size

% Matched filter's transfer function
H = sqrt(beta/pi)*exp(-1i*pi^2*fa.^2/beta);
H = repmat(H, rows, 1); % replicas of H on every row

% FFT on rows (azimuth) ---> (tau, fa) domain
image_tau_fa = fftshift(fft(ifftshift(raw_image, 2), [], 2), 2);

% azimuth compression
Y = image_tau_fa.*H;

% IFFT on rows (azimuth) ---> (tau, ta) domain
filtered_image = fftshift(ifft(ifftshift(Y, 2), [], 2), 2);

end