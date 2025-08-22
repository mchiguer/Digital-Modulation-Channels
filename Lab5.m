%% Ejercicio 1:

Nfft = 256;       
Ndata = 234;        
Nsym = 26;          
M = 64;       
Nb_bits = Ndata * Nsym * log2(M);   %Número total de bits transportados =234 * 26 * 6 bits
CP = 25; 
bits = randi([0 1], Nb_bits, 1); 
entiers = bit2int(reshape(bits, [], log2(M)), log2(M));
symbols = qammod(entiers, M, 'UnitAveragePower', true); 

X = zeros(Nfft, Nsym);  % Creamos una matriz de ceros de 256 * 26 = número de subportadoras * número de símbolos OFDM

symbol_matrix = reshape(symbols, Ndata, Nsym); 

X(2:Ndata/2+1,:) = symbol_matrix(1:Ndata/2,:);  % Rellenamos la parte superior de la matriz de ceros
X(Nfft - Ndata/2 + 1 : Nfft, :) = symbol_matrix(Ndata/2+1 : end, :); % Rellenamos la parte inferior de la matriz de ceros

X_OFDM = ifft(X, Nfft);

indice = [Nfft-CP+1:Nfft   1:Nfft]; %Añadimos las 25 filas con todas las columnas al inicio
X_con_PC = X_OFDM(indice,:);  
Vector_OFDM = reshape(X_con_PC, [], 1); 



%% Ejercicio 2: 

fs = 20e6;     % 20MHz
ts = 1 / fs;
fd = 5;        
SNR_dB = 25;

tau = [0, 0.3, 1.5, 3.4] * ts;   %Vector de retardos en segundos
gain_dB = [0, -1, -3, -2];       %Vector de ganancias

canal = comm.RayleighChannel('SampleRate', fs, 'ChannelFiltering', true, ...
    'PathDelays', tau, 'AveragePathGains', gain_dB, 'MaximumDopplerShift', fd); 

yofdm_rx = canal(Vector_OFDM); %Filtrado por el canal

yofdm_rx_awgn = awgn(yofdm_rx, SNR_dB, 'measured'); 



%% Ejercicio 3 : 

y_rx_cp = reshape(yofdm_rx_awgn, Nfft+CP, Nsym);

y_rx_noCP = y_rx_cp(CP+1:end, :); %Solo toma las líneas desde la 26 hasta el final

X_OFDM_trasmit_sin_CP = canal(X_OFDM(:,1));
X_OFDM_trasmit_sin_CP_Con_Ruido = awgn(X_OFDM_trasmit_sin_CP, SNR_dB,'measured');

figure;
subplot(2,1,1);
plot(abs(fft(X_OFDM_trasmit_sin_CP_Con_Ruido))); title('Espectro de símbolo OFDM transmitido (sin CP)');

subplot(2,1,2);
plot(abs(fft(y_rx_noCP(:,1)))); title('Espectro de símbolo OFDM recibido (después del canal + ruido)');

Y = fft(y_rx_noCP, Nfft, 1);  % FFT en cada columna



%% Ejercicio 4 : 


% Selección de símbolos recibidos después de FFT (Z) para pilotos
Z1 = Y(:,1);   
Z4 = Y(:,4);  
Z15 = Y(:,15); 

%Selección de símbolos transmitidos antes de IFFT (X) para los mismos índices
X1 = X(:,1);   
X4 = X(:,4);
X15 = X(:,15);

%Estimación porcentual del canal para cada controlador
H1 = Z1 ./ X1;
H4 = Z4 ./ X4;
H15 = Z15 ./ X15;


H_est = (H1 + H4 + H15) / 3; 

% Construcción de la matriz Hrep, que es una matriz con el vector de valor promedio
Hrep = repmat(H_est, 1, Nsym); 


% Ecualización elemento por elemento
Z_eq = Y ./ Hrep;  


%% Ejercicio 5 : 

% Supresión de portadora nula
Z_data = [Z_eq(2 : Ndata/2 + 1,:); Z_eq(Nfft - Ndata/2 + 1 : Nfft,:)];  

%conversión de Paralelo to serie
z_eq_vec = Z_data(:);

Z_raw = [Y(2 : Ndata/2 + 1, :); Y(Nfft - Ndata/2 + 1 : Nfft, :)]; 
z_raw_vec = Z_raw(:); 


%Visualización de constelación
figure;
subplot(1,2,1)
scatterplot(z_eq_vec);
title('Símbolos ecualizados');

subplot(1,2,2)
scatterplot(z_raw_vec);
title('Simbolos sin ecualización');


%% Ejercicio 6 : 

% Obtener bits
bits_rx = qamdemod(Z_received, M, ...
                   'OutputType', 'bit', ...
                   'UnitAveragePower', true); 

% Asegúrese de que los bits demodulados estén en formato de columna
bits_rx = bits_rx(:);

%Aplicar tasa de error entre bits transmitidos
errorCalc = comm.ErrorRate;          
resultats = errorCalc(bits, bits_rx);


% Mostrando resultados
fprintf('Número de bits con error: %d\n', resultats(2));
fprintf('BER estimado: %.5f\n', resultats(1));
fprintf('Número total de bits comparados: %d\n', resultats(3));



