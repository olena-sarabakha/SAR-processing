%%=========================================================================
%% HELP
%%[ENG]: function that plots the Point Spred Function (PSF) in slant range, with two
%%subplots dedicated to the visualization of the width at -4dB and the
%%level of the side lobes. Input: input_image, which must be a column of the data matrix;
%%swath_range, swath dimension in range in meters; ylim_lobo_laterale,
%%indicates the y limit to visualize the side lobes. Output: no outputs are
%%produced, but the function automatically opens a figure window and
%%displays the results revealed in intensity in dB and normalized to the
%%max value.


%%[ITA]: funzione che plotta la psf in slant range, con due subplot dedicati a larghezza a
%%-4dB e al livello dei lobi laterali. Input: input_image, immagine in
%%input, deve essere una colonna dell'immagine; swath_range, dimensione dello
%%swath in range, deve essere in metri; ylim_lobo_laterale, indica il ylim a cui
%%visualizzare i lobi laterale. Output: non sono presenti output, la
%%funzione apre automaticamente una funestra con l'immagine rivelata in
%%potenza in dB e normalizzata al suo valore massimo.
%%=========================================================================
function [] = psf_srange(input_image, swath_range, ylim_lobo_laterale)

[rows, ~] = size(input_image);
x = (-rows/2:rows/2-1)*swath_range/rows; % asse y (slant range) [m]

% generazione dell'immagine
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 2, [1, 2])
plot(x, 20*log10(abs(input_image)/max(abs(input_image))), 'k', 'LineWidth', 1.5)
text = sprintf('PSF in range \nPixel spacing=%f m', swath_range/rows);
title(text)
xlabel('Azimuth [m]')
ylabel('Intensità [dB]')
grid on
grid minor
subplot(2, 2, 3)
plot(x, 20*log10(abs(input_image)/max(abs(input_image))), 'k', 'LineWidth', 1.5)
ylim([-4, 0])
xlabel('Azimuth [m]')
ylabel('Intensità [dB]')
title('Apertura a -4dB')
grid on
grid minor
subplot(2, 2, 4)
plot(x, 20*log10(abs(input_image)/max(abs(input_image))), 'k', 'LineWidth', 1.5)
ylim([ylim_lobo_laterale, 0])
xlabel('Azimuth [m]')
ylabel('Intensità [dB]')
title('Lobi laterali')
grid on
grid minor

end