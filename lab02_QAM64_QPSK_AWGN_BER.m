%%2.1.
x_bits = randi([0 1], 5400, 1);

x_symbols = bit2int(x_bits,log2(64));

qam_modulated1 = qammod(x_symbols, 64); %valores enteros
qam_modulated2 = qammod(x_bits, 64, InputType = 'bit'); %bits

disp('Diferencia entre los dos valores :');
disp(norm(qam_modulated1 - qam_modulated2)); %para confirmar que coincidir

%%2.2. 

figure;
scatterplot(qam_modulated1(1:500), 1, 0, 'k*');
hold on;

scatterplot(qam_modulated2(1:500), 1, 0, 'mo');
hold off;

axis([-8 8 -8 8]); 
title('Constelación 64-QAM: primeros 500 símbolos');
grid on;

figure;
scatterplot(qam_modulated1(1:500), 100, 150, 'k*'); 
title('Constelación 64-QAM - Parámetros modificados');
grid on;


%%2.3. 
M = 64; 
x_bits = randi([0 1], 5400, 1);
x_symbols = bit2int(x_bits, log2(64)); 
y_qam_gray = qammod(x_symbols, 64, 'gray', 'PlotConstellation', true);
title('Constelación 64-QAM con codificación Gray');
figure;

y_qam_bin = qammod(x_symbols, 64, 'bin', 'PlotConstellation', true);
title('Constelación 64-QAM con mapeo de bits');
figure;

y_qam_bin_bits = qammod(x_bits, 64, 'InputType', 'bit', 'PlotConstellation', true);
title('Constelación 64-QAM: entrada binaria y mapeo binario');


%%3
%ej4 : 
bits_qpsk = randi([0 1], 10000, 1); 
qpsk_modulated = pskmod(bits_qpsk, 4, pi/4, 'InputType', 'bit');

canal_awgn = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)', ...
                              'EbNo', 18, 'BitsPerSymbol', log2(4));

P_signal = mean(abs(qpsk_modulated).^2);
canal_awgn.SignalPower = P_signal; 
qpsk_noisy = canal_awgn(qpsk_modulated);

figure;
subplot(1,2,1);
scatterplot(qpsk_modulated);
title('Constelación QPSK transmitida');

subplot(1,2,2);
scatterplot(qpsk_noisy);
title('Constelación QPSK con ruido (Eb/No = 18 dB)');


%ej5 : Cambia la constelacion
canal_awgn_reduit = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)', ...
                                     'EbNo', 6, 'BitsPerSymbol', log2(4));

qpsk_noisy_reduit = canal_awgn_reduit(qpsk_modulated);

figure;
subplot(1,2,1);
scatterplot(qpsk_noisy);
title('Constelación QPSK con ruido (Eb/No = 18 dB)');
subplot(1,2,2);
scatterplot(qpsk_noisy_reduit);
title('Constelación QPSK con ruido (Eb/No = 6 dB)');


%4
z = pskdemod(qpsk_noisy_reduit, 4, pi/4);

disp('Tipo de datos de z:');
disp(class(z)); 
disp('Primeros símbolos demodulados:');
disp(z(1:10)); 
disp('Primeros bits demodulados:');
disp(zbin(1:20));



%5
Error = comm.ErrorRate;

VectorErrores = Error(bits_qpsk, zbin);

BER = VectorErrores(1);

numero_erreurs = VectorErrores(2); 
disp(['BER : ', num2str(BER)]);
disp(['Número de bits erróneos: ', num2str(numero_erreurs)]);
