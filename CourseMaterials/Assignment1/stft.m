function [P, f, S] = stft(x, winparam, nfft, fs)
% short-time Fourier transform (STFT)
% // Input // %
% x:        the original data samples (Time Points x Channels)
% winparam: for an positive interger input, winparam is the window size (default window is hamming)
%           for a vector input, winparam is a window
% nfft:     number of fft points
% fs:       sampling rate
%
% // Output // %
% P:        squared magnitude of STFT (spectrogram)
% f:        evaluated frequency bins in STFT
% S:        complex time-frequency value of STFT

%% ======================================================= %%
% ELEC 6081 Biomedical Signals and Systems
% by Zhiguo Zhang, 09/2013
% ========================================================  %

fprintf('\nShort-time Fourier Transform: ')

%% Specify Parameters
if size(x,1)==1; x = x.'; end; % transpose data if the 1st dimension is 1
f = fs/2*linspace(0,1,nfft/2+1)';
N_Trials = size(x,2); % number of trials
N_T = size(x,1); % number of time samples
N_F = length(f); % number of frequency bins

fprintf('%d Time Points x %d Frequency Bins x %d Trial(s)\n',N_T,N_F,N_Trials);
S = single(zeros(N_F,N_T,N_Trials));
fprintf('Processing...     ')

%% Windowing and Padding
if length(winparam)==1 % a window size is specified
    if mod(winparam,2); h = winparam; %  window size (points)
    else h = winparam+1; %  window size (points); 
    end
    win = window('hamming',h); % window (one trial); default window type is hamming
else
    win = winparam;
    h = length(win);
end
W = repmat(win,1,N_Trials); % window (all trials)
U = win'*win;  % compensates for the power of the window

% Zero padding (default mode is "zero")
X = padarray(x,(h-1)/2);        % padding data
X = detrend(X);

%% STFT
for n=1:N_T
    fprintf('\b\b\b\b%3.0f%%',n/N_T*100)
    X_n = X(n+[0:h-1],:);
    X_n = detrend(X_n);
    S_n = fft(X_n.*W,nfft,1);
    S(:,n,:) = S_n(1:(nfft/2+1),:) / sqrt(U);
end
S2 = S.*conj(S);
P = S2/fs;
fprintf('  Done!\n')

end