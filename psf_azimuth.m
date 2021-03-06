%%=========================================================================
%% HELP
%%[ENG]: function that plots the Point Spred Function (PSF) in azimuth, with two
%%subplots dedicated to the visualization of the width at -4dB and the
%%level of the side lobes. Input: input_image, which must be a row of the data matrix;
%%Dsa, synthetic aperture of the system in meters; ylim_lobo_laterale,
%%indicates the y limit to visualize the side lobes. Output: no outputs are
%%produced, but the function automatically opens a figure window and
%%displays the results revealed in intensity in dB and normalized to the
%%max value.


%%[ITA]: funzione che plotta la psf, con due subplot dedicati a larghezza a
%%-4dB e al livello dei lobi laterali. Input: input_image, immagine in
%%input, deve essere una riga dell'immagine; Dsa, apertura sintetica del
%%sistema, deve essere in metri; ylim_lobo_laterale, indica il ylim a cui
%%visualizzare i lobi laterali. Output: non sono presenti output, la
%%funzione apre automaticamente una finestra con l'immagine rivelata in
%%potenza in dB e normalizzata al suo valore massimo.
%%=========================================================================

function [] = psf_azimuth(input_image, Dsa, ylim_lobo_laterale)

[~, columns] = size(input_image);
x = (-columns/2:columns/2-1)*Dsa/columns; % asse x (azimuth) [m]

% generazione dell'immagine
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2, 2, [1, 2])
plot(x, 20*log10(abs(input_image)/max(abs(input_image))), 'k', 'LineWidth', 1.5)
text = sprintf('PSF azimutale \nPixel spacing=%f m', Dsa/columns);
title(text)
xlabel('Azimuth [m]')
ylabel('Intensit� [dB]')
grid on
grid minor
subplot(2, 2, 3)
plot(x, 20*log10(abs(input_image)/max(abs(input_image))), 'k', 'LineWidth', 1.5)
ylim([-4, 0])
xlabel('Azimuth [m]')
ylabel('Intensit� [dB]')
title('Apertura a -4dB')
grid on
grid minor
subplot(2, 2, 4)
plot(x, 20*log10(abs(input_image)/max(abs(input_image))), 'k', 'LineWidth', 1.5)
ylim([ylim_lobo_laterale, 0])
xlabel('Azimuth [m]')
ylabel('Intensit� [dB]')
title('Lobi laterali')
grid on
grid minor

end