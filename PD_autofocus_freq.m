%%=========================================================================
%HELP: fuction that applies the Phase Difference Algorithnm (PD) in the FREQUENCY DOMAIN.
%%Input: input_image, in double time domain; R0, closest-approach
%%distance; Bd, Doppler bandwidth; ta, slow time vector; lambda, wavelength
%%of the signal. Output: estimated_beta, estimated value of the correct
%%focalization parameter; estimated_velocity, estimated value of the true velocity 
%%of the platform.

%tau = fast time; ta = slow time; fa = azimuth frequency.
%%=========================================================================

function [estimated_beta, estimated_velocity] = PD_autofocus_freq(input_image, ta, fa, Bd, lambda, R0)

[rows, columns] = size(input_image); % input matrix dimension

%% centering and windowing of the two looks in -B/4 e +B/4 and shifting into (tau, fa) domain
% look 1 in -B/4
look1 = input_image.*exp(-1i*2*pi*Bd/4*ta);
Look1 = fftshift(fft(ifftshift(look1, 2), [], 2), 2); % FFT on rows ---> (tau, fa) domain

% look 2 in +B/4
look2 = input_image.*exp(+1i*2*pi*Bd/4*ta);
Look2 = fftshift(fft(ifftshift(look2, 2), [], 2), 2); % FFT on rows ---> (tau, fa) domain

% windowing of the two looks
window = [zeros(rows, columns/4) ones(rows, columns/2) zeros(rows, columns/4)]; % window

Look1 = Look1.*window; % windowing of look 1
Look2 = Look2.*window; % windowing of look 2

%% multiplication of the two looks in (tau, fa) domain
Look1_conj = conj(Look1); % conjugate of look 1

DS = Look2.*Look1_conj;

%% estimate of beta in (tau, ta) domain
ds = ifftshift(ifft(fftshift(DS, 2), [], 2), 2); % IFFT on rows ---> (tau, ta)

% linear detection
ds_abs = abs(ds); 

% mean on columns
ds_mean = mean(ds_abs, 1); 

% peak detection
[~, max_pos] = max(ds_mean);

% estimated beta
estimated_beta = pi*Bd/(2*ta(max_pos));

% estimate of the true velocity
estimated_velocity = sqrt(lambda*R0*estimated_beta/(2*pi));

end















