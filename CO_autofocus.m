%%=========================================================================
%HELP: fucntion that applies the Contrast Optimization (CO) autofocus
%algorithm. Input: input_matrix, in double time domain; beta,
%nominal value of the focalizaion parameter; beta_ext, extention of beta
%values --> interval where to look for the correct value of beta; PRT, 
%pulse repetition time; lambda, wavelength of the signal; R0, closest-approach
%distance. Output: beta_est, estimated value of the focalization parameter;
%velocity_est, estimated true velocity; iter, number of iterations needed
%to achieve an accuracy such that the std of the phase is smaller than
%pi/4.

%tau = fast time; ta = slow time; fa = azimuth frequency.
%%=========================================================================

function [beta_est, velocity_est, iter]=CO_autofocus(input_matrix, beta, beta_ext, PRT, lambda, R0)

[rows, columns] = size(input_matrix); % input matrix dimension
ta = PRT*(-columns/2:columns/2-1); % slow time vector

beta_interval = [beta-beta_ext:1:beta+beta_ext]; % interval for beta values used in the attempts to find beta_est

n = 1; % initialization of the number of iterations
sigma_phi = inf; % initialization for sigma_phi
while sigma_phi > pi/4
    
    contrast = zeros(1, length(beta_interval)); % pre-allocation of contrast vector
    
    for b = 1:length(beta_interval)
        % dechirping 
         dechirp = exp(1i*beta_interval(b)*ta.^2); % reference for dechirping
         dechirp_matrix = repmat(dechirp, rows, 1); % matrix for dechirping
         matrix_dechirped = input_matrix.*dechirp_matrix; % dechirping of input matrix
         matrix_dechirped_fa = fftshift(fft(ifftshift(matrix_dechirped, 2), [], 2), 2); % FFT ---> (tau, fa)
         matrix_abs = abs(matrix_dechirped_fa); % abs value of the obtained matrix
         
         % contrast
         contrast(b) = std(matrix_abs, 1, 'all')/mean(matrix_abs, 'all');
%          contrast(b) = max(mean(matrix_abs, 1));
    end
    
    [~, max_index] = max(contrast); % find max value adn its position
    
    beta_aux(n) = beta_interval(max_index); % estimated beta in the current iteration
    dbeta_aux(n) = abs(beta_aux(n)-beta); % delta beta estimated in the current iteration
    
    % standard deviation of the phase
   if n == 1
    sigma_phi(n) = std(dbeta_aux(n)*ta.^2, 1, 'all');
       sigma_phi(n) = std((dbeta_aux(n)-dbeta_aux(n-1))*ta.^2, 1, 'all'); 
   end
   
   % stop criterion
   if sigma_phi(n) < pi/4
      break 
   end
  
   beta_interval = [beta_aux(n)-beta_ext/(n+1):1/(n+1):beta_aux(n)+beta_ext/(n+1)]; % update of beta interval where to look for the correct value
   n = n+1; % update of iteration index
end

beta_est = beta_aux(n); % final value of the estimated beta
velocity_est = sqrt(beta_est*lambda*R0/(2*pi)); % estimate of the true velocity
iter = n; % total number of iterations
end
       