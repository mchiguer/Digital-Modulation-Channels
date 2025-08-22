%% 2.
x = randi([0, 1], 1, 4000);
figure;
stem(x(1:10), 'filled');
title('First 10 bits of the numeric message');
xlabel('Index of the bit');
ylabel('Valor (0 ou 1)'); 
grid on;

%% 3.
baseband_signal = repelem(4 * x, 120);
t = (0:length(baseband_signal)-1) / 100e3;
figure;
plot(t(1:2000), baseband_signal(1:2000));
%Tb = Numero de muestras por bit / la frecuencia de muestra = 120 / 100e3
%R= 1/Tb = 100e3 / 120 

%%4.1
portadora = cos(2 * pi * 8e3 * t);
senal_modulada_ASK = baseband_signal .* portadora;
figure;
subplot(2,1,1);
plot(t(1:1000), baseband_signal(1:1000), 'b', 'LineWidth', 1.5); 
axis off; % Enlève les axes 
subplot(2,1,2);
plot(t(1:1000),senal_modulada_ASK(1:1000), 'r', 'LineWidth', 1.5);
axis off; % Enlève les axes
%T = N / fs = 1000/100e3;

%% 4.2
z=2*x-1; 
baseband_signal2= repelem(1 * z, 120);
figure;
plot(t(1:5*120), baseband_signal2(1:5*120), 'b', 'LineWidth', 1.5); 
portadora2 = cos(2 * pi * 8e3 * t);
senal_modulada_BPSK= baseband_signal2 .* portadora;
figure;
plot(t(1:1000),senal_modulada_BPSK(1:1000), 'b', 'LineWidth', 1.5); hold on;
figure;
plot(t(1:1000), baseband_signal2(1:1000), 'b', 'LineWidth', 1.5); 
%Los cambios de fase que detecta
phase_shifts = find(diff(z) ~= 0);
disp(phase_shifts);

%% 4.3 

window_size = 1024; 
noverlap = window_size / 2; % 50% de recouvrement
nfft = 4096;  
figure; 
pwelch(senal_modulada_BPSK, window_size, noverlap, nfft, fs);
title('DSP - Modulation BPSK');
figure;
pwelch(senal_modulada_ASK, window_size, noverlap, nfft, fs);
title('DSP - Modulation ASK');
%Las características más distintivas de ambos espectros:
% ASK: Picos principales alrededor de la frecuencia portadora (8 kHz) / Espectro más amplio porque la amplitud varía directamente con el mensaje.
% BPSK: Pico central alrededor de la frecuencia portadora / Espectro más estrecho que ASK, porque la potencia se concentra alrededor de la portadora.


%% 4.4 
symbols = randsample([-3, -1, 1, 3], 3000, true); %true c’est avec remise
% Séparation phase (I) et quadrature (Q) :
I_seq = symbols(1:2:end); % Indices impairs → phase 
Q_seq = symbols(2:2:end); % Indices pairs → quadrature

% Étendre chaque symbole :
I_signal = repelem(I_seq, 120); 
Q_signal = repelem(Q_seq, 120);
t = (0:length(I_signal)-1) /100e3;

samples_15ms = round(0.015 * fs);

figure;
subplot(2,1,1);
plot(t(1:samples_15ms), I_signal(1:samples_15ms), 'b', 'LineWidth', 1.5);
title('Composante en Phase (I)'); 
xlabel('Temps (s)'); 
ylabel('Amplitude');
grid on;
subplot(2,1,2);
plot(t(1:samples_15ms), Q_signal(1:samples_15ms), 'r', 'LineWidth', 1.5);
title('Composante en Quadrature (Q)'); 
xlabel('Temps (s)');
ylabel('Amplitude'); 
grid on;
portadora_cos = cos(2 * pi * 8e3* t);
portadora_sin = sin(2 * pi * 8e3* t);
QAM16_senal = I_signal .* portadora_cos- Q_signal .* portadora_sin;
figure;
plot(t(1:400), QAM16_signal(1:400), 'k', 'LineWidth', 1.5); 
title('QAM16_senal');

