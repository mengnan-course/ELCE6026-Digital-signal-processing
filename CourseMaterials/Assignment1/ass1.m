%% parameter

omiga_p = 0.6 * pi;
omiga_s = 0.7 * pi;
delta1 = 0.01;
delta2 = 0.001;

domiga = omiga_s - omiga_p;
delta = min(delta1,delta2);

A = -20*(log10(delta));

beta = [];
if A > 50
    beta = 0.1102 * (A-8.7);
elseif (A >= 21 && A <= 50)
    beta = 0.5842 * (A - 21)^0.4 + 0.07886 * (A - 21);
else
    beta = 0;
end

M = (A - 8)/(2.285*domiga);

M_pc = (-10*log10(delta1*delta2)-13)/(2.324*domiga);

%% (v)
n1 = ceil(M);
Wn = (omiga_p + omiga_s)/(2*pi);% 0 <= Wn <= 1
b1 = fir1(n1,Wn,kaiser(n1+1,beta));
fvtool(b1,1)

n2 = ceil(M_pc);
f = [0, omiga_p/pi, omiga_s/pi, 1];
m = [1 1 0 0];
b2 = firpm(n2,f,m);

fvtool(b1,1,b2,1);% b1 in blue; b2 in orange

%% (vii)
result1 = filter(b2,1,[x ; zeros(length(b2)+length(x)-1-length(x), 1)])
result2 = conv(b2,x)

all(abs(result1 - result2) < 0.001)

% plot result
plot(x);
title('The original signal');
figure;
plot(result1);
title ('The signal after filtering');

%% (viii)

% %% perform STFT
win_len=0.3; %% default setting for window length;
winsize =win_len*Fs; % actual window size for STFT;
nfft = 1024;   % # FFT points
[P, f] = stft(result1, winsize, nfft, Fs);

% %% display spectrogram
figure
imagesc(t,f,P)
colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
axis xy
grid on
set(gca,'ylim',[1 100]) % set the limits of frequency in the plot
title(['window size ' num2str(winsize)]);

