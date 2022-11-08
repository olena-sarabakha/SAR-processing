function [patch_filtered_image, dim_patch, beta] = patch_adapt_filter(raw_matrix, ps_range, da, lambda, altitude, L, va, fa, R0, beta_center)

%%=========================================================================
%%HELP: this function applies the matched filter in azimuth in a patched
%%fashion. The dimension of the patches is chosen in such a way to satisfy
%%the depth of focus criterion. Input: raw_matrix, I&Q matrix that has to
%%be focused; ps_range, pixel ranging in range; da, antenna's dimension in
%%azimuth; lambda, signal's wavelength; altitude, altitude of the SAR
%%platform; L, near range distance; va, velocity of the SAR platform; fa,
%%azimuth frequency vector; R0, closest-approach distance; beta_center,
%%focalization parameter relative to the center of swath. Output:
%%patch_filtered_image, filtered image; dim_patch, dimension of the
%%patches; beta, vector of focalization parameters for each patch.
%%=========================================================================

[rows, columns] = size(raw_matrix); % dimensions of the input matrix

distanza = zeros(rows, 1);

for cella_range=1:rows
    % distance of the cosidered range cell
    distanza(cella_range) = sqrt(altitude^2+(L+ps_range/2+(cella_range-1)*ps_range)^2); 
end

%% patch dimension definition
threshold = da^2/(2*lambda*R0); % threshold for beta variation allowed according to the depth of focus crioterion
d_beta = 0; % initialization value
i = rows/2+1;  

while d_beta/beta_center < threshold
    d_beta = abs((2*pi*va^2/(lambda*distanza(i)))-beta_center); % difference between the beta of the cosidered cell and beta of central cell (center of swath)
    i = i+1;
end

% patch dimension
dim_patch = 2*(i-2-rows/2)+1; % double of the obtained cells plus one (central cell)
% dim_patch = 3;

% number of patches
num_patch = floor(rows/dim_patch); 

%% matched filter filtering
beta = zeros(num_patch, 1); % preallocation of beta vector
patch_filtered_image = zeros(rows, columns); % preallocation of the output matrix
n = 1;

stop = num_patch*dim_patch;
for j=1:dim_patch:stop
        % beta of the patch
        beta(n) = 2*pi*va^2/(lambda*distanza(j+ceil(dim_patch/2)));
    
        % matched filtering
        [patch_filtrato] = matched_filter(raw_matrix(j:j+dim_patch, :), beta(n), fa);
    
        % filtered patch added to the output matrix
        patch_filtered_image(j:j+dim_patch, :) = patch_filtrato;

        n = n+1; % counter update
    
end

beta(n) = 2*pi*va^2/(lambda*distanza(rows)); % recording of beta values in a vector

% filtering of the last range cells, in case the patches don't cover the
% whole matrix
[patch_filtrato] = matched_filter(raw_matrix(j+dim_patch:end, :), beta(n), fa);
patch_filtered_image(j+dim_patch:end, :) = patch_filtrato;

end

    
