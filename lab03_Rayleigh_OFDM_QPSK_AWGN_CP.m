%%1.
%
fc = 2600e6;   % Fréquence porteuse (Hz)
ts = 130.2e-9; % Temps de muestreo (s)
c= 3e8;
v1 = 300 *(1000/3600); % vitesse du canal en m/s
v2= 5 * (1000/3600); % vitesse du canal en m/s

% frecuencias Doppler:
fd1 =   v1 * fc / c  ; % Convertir en m/s la vitesse du canal 
fd2 = v2 * fc / c ; 

%Diseno del canales : 

% Channel c1 (2 taps)
channel_c1 = comm.RayleighChannel( ...
    'SampleRate', 1/ts, ...
    'PathDelays', [0, 0.3 * ts], ...
    'AveragePathGains', [0, -5], ...
    'MaximumDopplerShift', fd1);

% Channel c2 (5 taps)
channel_c2 = comm.RayleighChannel( ...
    'SampleRate', 1/ts, ...
    'PathDelays', [0, 0.9*ts, 1.6*ts, 1.9*ts, 2.3*ts], ...
    'AveragePathGains', [0, -1, -3, -2, -4], ...
    'MaximumDopplerShift', fd2);

c1.ChannelFiltering = true; 
c2.ChannelFiltering = true; 


% Genera una secuencia de 1 000 bits.
num_bits = 1000;
bits = randi([0 1], 1, num_bits);

% agrupar los bits en pares
bit_pairs = reshape(bits, 2, []); % Matriz de 2 filas y (num_bits/2) columnas

symbols = bit2int(bit_pairs', 2);  % convertir en vecteur colonne

modulated_qpsk = pskmod(symbols, 4, pi/4); 

% Filtrar la señal modulada con ambos canales
received_c1 = channel_c1(modulated_qpsk);
received_c2 = channel_c2(modulated_qpsk);


% Representar las constelaciones resultantes
figure;
scatterplot(received_c1);
title('Constelación recibida - Canal c1 (300 km/h)');

figure;
scatterplot(received_c2);
title('Constelación recibida - Canal c2 (5 km/h)');

% Comentario 
%El canal c1 es selectivo en frecuencia y rápido,
% mientras que el canal c2 es casi plano en frecuencia.

%%2. 
N_IFFT = 256;
N_null = 46;
N_data = N_IFFT - 2*N_null - 2;

% Un vector que contenga Ndata símbolos QPSK 
qpsk_data_used = modulated_qpsk(1:N_data); 

% Diseña un vector Z
Z = zeros(N_IFFT, 1);
Z(2:N_data/2+1) = qpsk_data_used(1:N_data/2);                   % frecuencias positivas
Z(end - N_data/2 + 1:end) = qpsk_data_used(N_data/2+1:end);     % frecuencias negativas

% Genera el símbolo OFDM z[n] mediante la transformada iFFT 
z_t = ifft(Z);


%visualiza el espectro del simbolo OFDM 
[PSD_tx, f_tx] = pwelch(z_t, [], [], [], 1);
figure; 
plot(f_tx, 10*log10(PSD_tx)); 
grid on;

% Transmite el símbolo OFDM por el canal c2 
z_rayleigh = channel_c2(z_t); 

%visualiza el espectro de la señal recibida
[PSD_rayleigh, f_rayleigh] = pwelch(z_rayleigh, [], [], [], 1);
figure; 
plot(f_rayleigh, 10*log10(PSD_rayleigh)); 
grid on;

% Añade a la señal recibida un ruido blanco gaussiano 
canal_awgn = comm.AWGNChannel( ...
    'NoiseMethod', 'Signal to noise ratio (Eb/No)', ...
    'EbNo', 15, ...
    'BitsPerSymbol', 2);

z_awgn = canal_awgn(z_rayleigh);
[PSD_awgn, f_awgn] = pwelch(z_awgn, [], [], [], 1);
figure; 
plot(f_awgn, 10*log10(PSD_awgn)); 
grid on;

%El ruido AWGN afecta la señal en todas las frecuencias. 
% A diferencia del canal de Rayleigh, que actúa de forma selectiva en frecuencia.


% Un prefijo cíclico de longitud NCP = 16
N_CP = 16; 

% Añadir un prefijo cíclico de longitud NCP = 16 al símbolo a transmitir
z_cp = [z_t(end-N_CP+1:end); z_t];

[PSD_nocp, f1] = pwelch(z_t, [], [], [], 1);
[PSD_cp, f2]   = pwelch(z_cp, [], [], [], 1);

figure;
plot(f1, 10*log10(PSD_nocp)); 
hold on;
plot(f2, 10*log10(PSD_cp)); 
grid on;

%Se observa que los espectros con y sin prefijo son prácticamente idénticos,
% porque el prefijo cíclico actúa sólo en el dominio del tiempo para combatir la interferencia entre símbolos (ISI).
