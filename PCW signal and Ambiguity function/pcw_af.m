%%%%%
%%% Calculate the auto-ambiguity function of single-period and four-period PCW signals
%%%
clc
clear
close all
fs = 50e3;     % sampling frequency 50kHz
T = 0.1;       % duration
B = 1/T;       % bandwidth
p_wid = 0.02;  % single-period pulse width
pcw4_t = 0.005;% four-period pulse width
t = 0:1/fs:T-1/fs;
t1 = 0:1/fs:p_wid-1/fs;
t_pcw4 = 0:1/fs:pcw4_t-1/fs;
prd = 0.025;         % pulse repetition interval
prd_sam = 0.025*fs;  % Pulse repetition interval sampling points
A_pcw = 1;      % amplitude
f_pcw = 7e3;    % center frequency 7kHz
I = 4;          % pulse number
c = 1500;       % speed of sound
% R = [1025,1040,1053,1060];
% v = [5,-8,8,-10];
% rou = [0.8,1,0.7,0.6];
% v_dopl = 2*f_pcw/c.*v;
y_pcw1 = A_pcw*cos(2*pi*f_pcw*t1);
pcw4_1 = A_pcw*cos(2*pi*f_pcw*t_pcw4);
y_pcw = [y_pcw1,zeros(1,length(t)-length(t1))];   % generate a single-period PCW  signal
pcw4 = zeros(1,length(t));                        % a four-period PCW signal is generated as follows
for i = 1:I
    pcw4((i-1)*prd_sam+1:(i-1)*prd_sam+length(pcw4_1)) = pcw4_1;
end

N = 2^nextpow2(fs*T);         % DFT points
Y = fft(y_pcw,N)/N*2;         % single-period PCW signal spectrum
Y_pcw4 = fft(pcw4,N)/N*2;     % four-period PCW signal spectrum
f = fs/N*(0:1:N-1);           % Frequency axis
P = abs(Y);                   % single-period PCW signal magnitude spectrum
P_pcw4 = abs(Y_pcw4);         % four-period PCW signal magnitude spectrum
maxDelay = 0.1;               % maximum delay of calculation
maxDoppler = 1000;            % calculated maximum Doppler frequency shift
[af_pcw,delay,doppler] = computeAmbiguityFunction(y_pcw,fs,maxDoppler,maxDelay);
[af_pcw4,delay_4,doppler_4] = computeAmbiguityFunction(pcw4,fs,maxDoppler,maxDelay);
af_pcw_d = [fliplr(af_pcw'),af_pcw'];
delay_d = [-fliplr(delay),delay];
a_nom_max = 0;                % find the maximum value below
for i = 1:length(af_pcw_d(:,1))
    a_nom_max = max(a_nom_max,max(af_pcw_d(i,:)));
end
af_pcw_d_nom = af_pcw_d./a_nom_max;  % normalize
af_pcw4_d = [fliplr(af_pcw4'),af_pcw4'];
delay_4_d = [-fliplr(delay_4),delay_4];
a_nom_max = 0;                % find the maximum value below
for i = 1:length(af_pcw4_d(:,1))
    a_nom_max = max(a_nom_max,max(af_pcw4_d(i,:)));
end
af_pcw4_d_nom = af_pcw4_d./a_nom_max;% normalize
%%
figure
subplot(211)
plot(t,y_pcw);
xlabel('t/s'); ylabel('magnitude/v'); title('time domain waveform');
subplot(212)
plot(f(1:N/2+1),P(1:N/2+1));
xlabel('frequency(Hz)');  ylabel('magnitude');  title('frequency spectrum');
figure
subplot(211)
plot(t,pcw4);
xlabel('t/s'); ylabel('magnitude/v'); title('time domain waveform');
subplot(212)
plot(f(1:N/2+1),P_pcw4(1:N/2+1));
xlabel('frequency(Hz)');  ylabel('magnitude');  title('frequency spectrum');
%%
figure
mesh(delay_d,doppler,af_pcw_d_nom);
colormap jet
colorbar
xlabel('time delay (s)');
ylabel('Doppler frequency shift (Hz)');
title('single PCW Auto-AF');
figure
mesh(delay_4_d,doppler_4,af_pcw4_d_nom);
colormap jet
colorbar
xlabel('时间延迟 (s)');
ylabel('多普勒频移 (Hz)');
title('fc=7k,num=4,t=5ms,T=25ms four PCW-Auto-AF');
