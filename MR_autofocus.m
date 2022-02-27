
%%=========================================================================
%%HELP: fuction that applies the Multilook Registration Algorithm to
%%autofocus the image. Input: input_image, in double time domain;
%%Tsa, synthetic aperture time; beta, nominal value of the ficalization
%%parameter; lambda, wavelength of the signal; R0, closest-approach
%%distance; PRF, pulse repetition frequency; da, antenna dimension in
%%azimuth. Output: betaMR, focalization parameter estimated by the
%%algorithm; Va_MR, estimated platform velocity; Tsa_MR, estimated
%%synthetic aperture time; iter, number of iterations needed to achieve an
%%accuracy such that the std of the phase is smaller than pi/4.

%tau = fast time; ta = slow time; fa = azimuth frequency.
%%=========================================================================

function [betaMR, Va_MR, Tsa_MR, iter]=MR_autofocus(input_image, ta, Tsa, fa, beta, lambda, R0, PRF, da)

[rows, columns] = size(input_image);
columns_new = 4*(2^nextpow2(columns)); % new number of colums for zero padding

ta_new = (1/PRF)*[-columns_new/2:columns_new/2-1]; % new slow time vector
fa_new = (PRF/columns_new)*[-columns_new/2:columns_new/2-1]; % new azimuth frequency vector

sigma_phi = inf; % initialization value for the cycle

i = 1; % counter initialization

% zero padding in azimuth domain in (tau, ta)
input_image = [zeros(rows, (columns_new-columns)/2) input_image zeros(rows, (columns_new-columns)/2)]; 

% nominal dechirping
b = exp(1i*beta*ta_new.^2); % reference chirp for dechirping
% b=repmat(b, R, 1); % matrice per dechirping
image_dechirped = input_image.*b; % nominal dechirping 

    while sigma_phi>pi/4
    
        % selection of the two looks
    look1 = image_dechirped;
    % figure, plot(abs(look1(26, :)))
    % hold on
    look2 = image_dechirped;
    % plot(abs(look2(26, :)))
    % hold off
    % legend('look1', 'look2')

        % FFT on rows ---> (tau, fa)
    Look1 = fftshift(fft(ifftshift(look1, 2), [], 2), 2);
    Look2 = fftshift(fft(ifftshift(look2, 2), [], 2), 2);

        % translation of the two looks
    Look1 = Look1.*exp(-1i*2*pi*Tsa/4*fa_new); % look translated in -Tsa/4
    % figure, plot(abs(Look1(26, :)))
    % hold on

    Look2 = Look2.*exp(+1i*2*pi*Tsa/4*fa_new); % look translated in +Tsa/4
    % plot(abs(Look2(26, :)))
    % hold off
    % legend('Look1', 'Look2')

        % windowing
    window = [zeros(rows, columns_new/2-columns/4) ones(rows, columns/2) zeros(rows, columns_new/2-columns/4)];
    % figure, plot(window(26, :))
    % hold on

        % IFFT on rows ---> (tau, ta)
    look1 = fftshift(ifft(ifftshift(Look1, 2), [], 2), 2);
    look2 = fftshift(ifft(ifftshift(Look2, 2), [], 2), 2);

        % windowing in (tau, ta) domain
    look1 = look1.*window;
    look2 = look2.*window;

        % FFT ---> (tau, fa)
    Look1 = fftshift(fft(ifftshift(look1, 2), [], 2), 2);
    Look2 = fftshift(fft(ifftshift(look2, 2), [], 2), 2);

        % linear detection on rows
    % for r=1:rows  
    %     abs_Look1(r, :) = abs(Look1(r,:));
    %     abs_Look2(r, :) = abs(Look2(r,:));
    % end
    abs_Look1 = abs(Look1);
    abs_Look2 = abs(Look2); 

        % mean value in order to obtain vectors
    mean1 = mean(abs_Look1, 1);
    mean2 = mean(abs_Look2, 1);

    % figure, plot(mean1)
    % hold on
    % plot(mean2)
    % legend('mean1', 'mean2')

        % cross-correlation to find peak
    corr = xcorr(mean1, mean2);
    % [r_corr, c_corr]=size(corr)
    % figure, plot(corr)
    f_corr = (PRF/columns_new)*[-columns_new:columns_new-2]; % frequency vector for correlation
    % [r_f, c_f]=size(f_corr)

    [~, ind_picco] = max(corr); % find max of correlation
    dbeta(i) = f_corr(ind_picco)*2*pi/Tsa; % estimate of delta_beta from the position of the peak in trequency vector
    beta_stim(i) = beta+sum(dbeta); % sum of delta_beta estimated in the current iteration to the previous estimates

    sigma_phi = rms(dbeta(i)*ta_new.^2);

    if sigma_phi<pi/4
        break
    end

    correz = exp(1i*beta_stim(i)*ta_new.^2);
    correz = repmat(correz, rows, 1);

    % image=fftshift(fft(ifftshift(image, 2), [], 2), 2);
    image_dechirped = input_image.*correz;
    % image=ifftshift(ifft(fftshift(image, 2), [], 2), 2);


    i = i+1;

    end

betaMR = beta_stim(i); % final estimated beta value   
Va_MR = sqrt(lambda*R0*betaMR/(2*pi)); % esrtimate of the true velocity
Tsa_MR = lambda*R0/(Va_MR*da);
iter = i; % total number of iterations
end