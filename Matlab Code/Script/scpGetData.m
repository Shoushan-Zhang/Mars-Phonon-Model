clear all
close all
clc

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig));
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'));

% Build a lookup table for 8 impacts
matchtable(1).ID="S0533a";matchtable(1).Date=datetime(2020,05,27,13,43,03);
matchtable(2).ID="S0793a";matchtable(2).Date=datetime(2021,02,18,19,31,23);
matchtable(3).ID="S0981c";matchtable(3).Date=datetime(2021,08,31,03,59,01);
matchtable(4).ID="S0986c";matchtable(4).Date=datetime(2021,09,05,05,18,58);
matchtable(5).ID="S1000a";matchtable(5).Date=datetime(2021,09,18,17,46,20);
matchtable(6).ID="S1034a";matchtable(6).Date=datetime(2021,10,23,18,26,30);
matchtable(7).ID="S1094b";matchtable(7).Date=datetime(2021,12,24,22,38,02);
matchtable(8).ID="S1160a";matchtable(8).Date=datetime(2022,03,02,06,51,40);


load sol_1034.mat % load seismic data for sol_XXX

% acceleration vertical Z: e.vbbzne(1).a
% acceleration north/south N: e.vbbzne(2).a
% acceleration east/west E: e.vbbzne(3).a

srSeis = e.vbbsrmax;     % sampling rate
% ------------- %
%
figure(7)
clf
%subplot(2,1,2)
T_interval = 10;              % Create frame with 10 seconds
sampleNumber = T_interval*srSeis;           % Each frame contains 10*20=200 samples
w  = hann(floor(sampleNumber)); % hanning window
N_window = 1.05; % overlap window for spectrogram calculation = 1/N_window
% Remove long-term trend, high-pass, get spectrogram for total power
% between all three components of SEIS ZNE
[~, seis.f, seis.t, seis.Z] = (spectrogram(SPremovelow(e.vbbzne(1).a,srSeis,1/T_interval), w, ...
    floor(sampleNumber/N_window), ...
    floor(sampleNumber/N_window), ...
    srSeis,'yaxis'));
[~, seis.f, seis.t, seis.N] = (spectrogram(SPremovelow(e.vbbzne(2).a,srSeis,1/T_interval), w, ...
    floor(sampleNumber/N_window), ...
    floor(sampleNumber/N_window), ...
    srSeis,'yaxis'));
[~, seis.f, seis.t, seis.E] = (spectrogram(SPremovelow(e.vbbzne(3).a,srSeis,1/T_interval), w, ...
    floor(sampleNumber/N_window), ...
    floor(sampleNumber/N_window), ...
    srSeis,'yaxis'));
seis.t = e.time(1) + seconds(seis.t);     % Sum up the power of all directions
seis.totalPower = seis.Z + seis.N + seis.E; %can use total power or single component, E, N or Z
%subplot(2,1,1)
surf(seis.t, seis.f, log10(sqrt(seis.totalPower)), 'EdgeColor', 'none')
shading(gca, 'flat');
set(gca,'FontSize',18)
grid off
colormap(turbo)
view(0,90)
caxis([-9.5 -6.5])
set(gca,'YScale','lin'); %set axis to linear to make visualisation faster
h = colorbar
ylabel(h, 'Acceleration (log_{10}(m/s^2/Hz^{1/2}))')
ylabel('Frequency (Hz)')
%% Get narrow-banded envelopes, for example the lander resonance around 4 Hz
% ------------- %
%subplot(2,1,2)
fmin = 3.7;
fmax = 4.3;
idxSeis(1) = find(seis.f > fmin,1,'first');
idxSeis(2) = find(seis.f < fmax,1,'last');
dfSeis = seis.f(3) - seis.f(2); % frequency resolution of spectrogram
hiFreq =  (sqrt(sum(seis.totalPower(idxSeis(1):idxSeis(2),:))*dfSeis));

% plot power of pressure data between 0.1 to 4 Hz;
figure(234)
clf
ylabel('Acceleration (log_{10}(m/s^2 Hz^{-1/2}))')
set(gca,'FontSize',18)
fmin2 = 0.02;
fmax2 = 1;
idxSeis(1) = find(seis.f > fmin2,1,'first');
idxSeis(2) = find(seis.f < fmax2,1,'last');
lowFreq = (sqrt(sum(seis.totalPower(idxSeis(1):idxSeis(2),:))*dfSeis)); % get low frequency for example 0.02 - 1 Hz
plot(seis.t,hiFreq,'LineWidth',2)
ylabel('Acceleration (m/s^2)')
set(gca,'FontSize',18)
title(['Envelope of seismic total power between ', num2str(fmin),' - ',num2str(fmax),' Hz'])



%%  S1034a

figure(234)
xlim([matchtable(6).Date-minutes(20),matchtable(6).Date+minutes(30)])

figure(7)
xlim([matchtable(6).Date-minutes(20),matchtable(6).Date+minutes(30)])

% Normalise the data



