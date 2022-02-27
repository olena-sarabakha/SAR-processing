function [output_image]=RCM_correction(input_image, R0, lambda, va, fa, fr)

%%=========================================================================
%%[ENG]
%%HELP: function that corrects the Range Cell Migration (RCM). The
%%correcting term deltaR is estimated from the acquisition geometry and
%%some parameters of the system. Input: input_image, image in double time
%%domain; R0, closest-approach distance; lambda, wavelength of the
%%signal; va, the system's velocity; fa, frequency vector in azimuth; fr,
%%frequency vector in range. Output: output_image, image with RCM
%%correction.

%%tau = fast time; ta = slow time; fr = range frequency; fa = azimuth
%%frequency.

%%[ITA]
%%HELP: funzione che corregge la Range Cell Migration (RCM). Il termine
%%correttivo deltaR viene stimato in base alla geometria di acquisizione e
%%ad alcuni parametri che caratterizzano il sistema. Input: input_image,
%%immagine da correggere; R0, distanza al closest approach; lambda,
%%lunghezza d'onda utilizzata; va, velocità della piattaforma; fa,
%%vettore della frequenza in azimuth; fr, vettore della frequenza in range.
%%Output: output_image, immagine con RCM compensata.
%%=========================================================================

c = 3e8; % speed of light [m/s]

% FFT to double frequency domain (fr, fa)
image_fr_ta = fftshift(fft(ifftshift(input_image, 1), [], 1), 1); % --> (fr, ta) domain   
image_fr_fa = fftshift(fft(ifftshift(image_fr_ta, 2), [], 2), 2); % --> (fr, fa) domain

% correction term deltaR and correction matrix
deltaR = R0*lambda^2/(8*va^2)*fa.^2; 
corr_matrix = exp(1i*4*pi*fr'*deltaR/c); 

% RCM correction
output_image = image_fr_fa.*corr_matrix;

% IFFT to go back to double time domain
output_image = fftshift(ifft(ifftshift(output_image, 1), [], 1), 1);
output_image = fftshift(ifft(ifftshift(output_image, 2), [], 2), 2);

end