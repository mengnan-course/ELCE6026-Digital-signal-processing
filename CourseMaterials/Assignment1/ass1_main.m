clear all; 
close all;

%% ======================================================= %%
% ELEC 2204 Digital Signal Processing
% Assignment 1 (Ver: Oct. 22, 2015)
% Prof. S. C. Chan, Dr. H. C. Wu 
% Special thanks to Dr. Z. G. Zhang for providing the signals and codes for analysis.


%% Load the synthetic signal
load Signal.mat
Fs=200; %% Sampling Frequency
Ts=1/Fs; %% Sampling Period
t=0:Ts:2;
t=t';
%% Plot the signal
figure
plot(t,x)
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-1.5 1.5])

%% perform STFT
win_len=0.3; %% default setting for window length;
winsize =win_len*Fs; % actual window size for STFT;
nfft = 1024;   % # FFT points
[P, f] = stft(x, winsize, nfft, Fs);

%% display spectrogram
figure
imagesc(t,f,P)
colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
axis xy
grid on
set(gca,'ylim',[1 100]) % set the limits of frequency in the plot
title(['window size ' num2str(winsize)]);



