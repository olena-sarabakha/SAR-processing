function [beta_estim, dVa]=PD_autofocus_time(image, fa, R0, lambda, Tsa, ta, beta, PRF)

%%=========================================================================
%%HELP: fuction that applies the Phase Difference Algorithnm (PD) in the TIME DOMAIN.
%%Input: image, image that has to be focused; R0, closest-approach
%%distance; beta, nominal value of the focalization parameter; PRF, pulse
%%repetition frequency. Output: beta_estim, estimated value of the correct
%%focalization parameter; dVa, estimated value of the true velocity of the
%%platform.

%%!!! This algorithm is NOT SUITABLE in case multiple scatterers are
%%present in the same range cell.


%%tau = fast time; ta = slow time; fr = range frequency; fa = azimuth frequency.

%%=========================================================================


[rows, columns] = size(image);

columns_new = 4*(2^nextpow2(columns)); % new number of columns for the padding (4 is an arbitrary multiplier)
% zero padding (tau, ta) domain
image = [zeros(rows, (columns_new-columns)/2) image zeros(rows, (columns_new-columns)/2)];

% new vectors for slow time and azimuth frequency
ta_new = (1/PRF)*[-columns_new/2:columns_new/2-1];
fa_new = (PRF/columns_new)*[-columns_new/2:columns_new/2-1];

b = exp(1i*beta*(ta_new.^2)); % exp term for a nominal correction
image = image.*b; % dechirping (battimento)

% division into two look

Look1 = fftshift(fft(ifftshift(image, 2), [], 2), 2); % look 1 in (tau, fa) domain
Look1 = Look1.*exp(-1i*2*pi*Tsa/4*fa_new); % centering of look 1
look1 = ifftshift(ifft(fftshift(Look1, 2), [], 2), 2); % look 1 centered in (tau, ta) domain

Look2 = fftshift(fft(ifftshift(image, 2), [], 2), 2); % look 2 in (tau, fa) domain
Look2 = Look2.*exp(1i*2*pi*Tsa/4*fa_new); % centering look 2
look2 = ifftshift(ifft(fftshift(Look2, 2), [], 2), 2); % look 2 centered in (tau, ta) domain


% windowing of the two looks in order to select the useful portion
window = [zeros(rows, columns_new/2-columns/4) ones(rows, columns/2) zeros(rows, columns_new/2-columns/4)];
look1 = look1.*window;
look2 = look2.*window;

% conj of look 2
    for r=1:rows
        look2_conj(r, :) = conj(look2(r, :));
    end

% multiplication of the two looks
ds = look1.*look2_conj;


% FFT in azimuth --> (tau, fa) domain
DS = fftshift(fft(ifftshift(ds, 2), [], 2), 2);


% linear detection
    for r=1:rows
        abs_DS(r, :) = abs(DS(r, :));
    end

% average on columns
av_DS = mean(abs_DS, 1);

% peak value detection
[~, pos_max] = max(av_DS);

% estimate of delta beta
dbeta = fa_new(pos_max)*2*pi/Tsa;
beta_estim = beta+dbeta;

dVa = sqrt(lambda*R0*beta_estim/(2*pi)); % estimate of correct velocity

end