%%1.

%1
p = [0.6, 0.35, -0.5];  % Canal

% la respuesta impulsional
figure;
stem(p);
title('Respuesta impulsional del canal');

% la respuesta en frecuencia
fvtool(p, 1);

%Si, se trata de un canal selectivo en frecuencia 
%porque atenúa de forma diferente las distintas componentes frecuenciales del señal.


%2 
Q = 2;
Pc = toeplitz([p zeros(1, Q)], [p(1) zeros(1, Q)]);
disp('Matriz Pc:');
disp(Pc);
%(Q+1)×(K+Q+1)=(2+1)×(2+2+1) = 3 × 5 son las dimensiones de la matriz

%3 
%Dimensiones de c son :  (K+Q+1) × 1  = (2+2+1) × 1 = 5 × 1
d = 1; 
K=2;

L = K + Q + 1;        %longitud del canal global
c = zeros(L, 1);      % inicializar vector c
c(d+1) = 1;           % ubicar 1 en la posición d+1
disp('Vector de canal ideal c:');
disp(c);

%4
wZF = (Pc' * Pc) \ (Pc' * c); % or : inv(Pc.'*Pc)*Pc.'*c; 
​disp('los coeficientes del filtro igualador :');
disp(wZF);


%5 
canal_global = conv(p, wZF.');  %Convolucion del canal et del igualador ZF
disp('Respuesta impulsional del canal global (p * wZF) :');
disp(canal_global);

%6 

p = [0.6, 0.35, -0.5];          
wZF = [0.5105; 0.7865; 0.1715];   
canal_global = conv(p, wZF.'); 


[Hc, w] = freqz(p, 1, 512);  
[Heq, ~] = freqz(wZF, 1, 512);          
[Htot, ~] = freqz(canal_global, 1, 512);  


figure;
axis([-1 1 -10 10])
plot(w/pi, 20*log10(abs(Hc)), 'b', 'LineWidth', 1.5); hold on;
plot(w/pi, 20*log10(abs(Heq)), 'r', 'LineWidth', 1.5);
plot(w/pi, 20*log10(abs(Htot)), 'g', 'LineWidth', 1.5);

grid on;
xlabel('Frecuencia normalizada (\times\pi rad/muestra)');
ylabel('Magnitude (dB)');
legend('Canal', 'Igualador ZF', 'Canal global');
title('Magnitude Response');


%%2

%parametros
p = [0.6, 0.35, -0.5];   
Q = 30;                 
K = 2;    %length(p)-1
d = 20;        

% 1
Pc = toeplitz([p, zeros(1, Q)], [p(1), zeros(1, Q)]);

%  2
c = zeros(K + Q + 1, 1);  
c(d + 1) = 1;             
disp(['tamano de c : ', num2str(length(c))]);

% 3
wZF = (Pc' * Pc) \ (Pc' * c);   % or : inv(Pc.'*Pc)*Pc.'*c; 
disp(['Número de coeficientes de filtro: ', num2str(length(wZF))]);

% 4
canal_global = conv(p, wZF.');

% 5
[Hc, w] = freqz(p, 1, 512);            
[Heq, ~] = freqz(wZF, 1, 512);          
[Htot, ~] = freqz(canal_global, 1, 512);  

figure;
plot(w/pi, 20*log10(abs(Hc)), 'b', 'LineWidth', 1.5); hold on;
plot(w/pi, 20*log10(abs(Heq)), 'g', 'LineWidth', 1.5);
plot(w/pi, 20*log10(abs(Htot)), 'r', 'LineWidth', 1.5); hold off;
grid on;
xlabel('Frecuencia normalizada (\times\pi rad/muestra)');
ylabel('Magnitud (dB)');
legend('Canal', 'Igualador ZF', 'Canal Igualado');
title(['Respuesta en frecuencia con Q = ' num2str(Q) ' y d = ' num2str(d)]);
axis([-1 1 -20 20]);


%%3 

% Disena una senal digital compuesta por 500 simbolos modulados en QPSK
NSym = 500;   
M = 4;% Numero de simbolos
bitsPerSym = log2(M);         % Numero de bits per simbolo
Nbits = NSym * bitsPerSym;      % Total_bits = 1000
                        

bits = randi([0 1], Nbits, 1);

symbols = pskmod(bit2int(bits, bitsPerSym), M, pi/4);

% 1
p = [0.6, 0.35, -0.5];  
x_channel = filter(p, 1, symbols);

% 2
SNR_dB = 15;         
x_noisy = awgn(x_channel, SNR_dB, 'measured');

%3
figure;
subplot(1,2,1);
scatterplot(x_channel);
title('Constelacion con ISI (canal)');
subplot(1,2,2);
scatterplot(x_noisy);
title('Constelacion con ISI + ruido');

% 4
Q = 30;
d = 20; 

Pc = toeplitz([p, zeros(1,Q)], [p(1), zeros(1,Q)]);

c = zeros(K+Q+1,1); 
c(d+1) = 1;

wZF = (Pc' * Pc) \ (Pc' * c);
x_equalized = filter(wZF, 1, x_noisy); 

% 5
figure;
subplot(1,2,1);
scatterplot(x_noisy);
title('Senal con bruit y ISI)');
subplot(1,2,2);
scatterplot(x_equalized);
title('Senal tras el ZF');

%Antes de ecualizar, los puntos están muy dispersos por el canal y el ruido.
%Después de aplicar el ZF, los puntos se juntan mejor en las 4 esquinas del QPSK.
%El ecualizador ayuda a corregir la distorsión del canal.



%%4 

%Parametros
p = [0.6, 0.35, -0.5];  
K = 2;        
Q = 30;               
d = 19;     

Pc = toeplitz([p, zeros(1, Q)], [p(1), zeros(1, Q)]);

c = zeros(K + Q + 1, 1); 
c(d + 1) = 1;        

% ZF
wZF = (Pc' * Pc) \ (Pc' * c);

[H_ZF, w] = freqz(wZF, 1, 512);

% Valores de SNR
SNR_dBs = [25, 15, 5];

colores = ['r', 'g', 'm'];  

% 
figure;
plot(w/pi, 20*log10(abs(H_ZF)), 'b', 'LineWidth', 1.5); hold on;
for i = 1:length(SNR_dBs)
    SNR_lin = 10^(SNR_dBs(i) / 10);
    lambda = 1 / SNR_lin;
    %MMSE
    wMMSE = (Pc' * Pc + lambda * eye(Q + 1)) \ (Pc' * c);
    [H_MMSE, ~] = freqz(wMMSE, 1, 512);
    plot(w/pi, 20*log10(abs(H_MMSE)), colores(i), 'LineWidth', 1.5);
end
xlabel('Frecuencianormalizada (\times\pi rad/muestra)');
ylabel('dB');
legend('ZF', 'MMSE 25dB', 'MMSE 15dB', 'MMSE 5dB');
title('Comparación ZF y MMSE (solo filtros)');
grid on;

% canal*igual
figure;
for i = 1:length(SNR_dBs)
    SNR_lin = 10^(SNR_dBs(i) / 10);
    lambda = 1 / SNR_lin;

    %MMSE
    wMMSE = (Pc' * Pc + lambda * eye(Q + 1)) \ (Pc' * c);
    canal_global = conv(p, wMMSE.');
    [Htot, w] = freqz(canal_global, 1, 512);

    % CanalGlobal
    plot(w/pi, 20*log10(abs(Htot)), colores(i), 'LineWidth', 1.5); hold on;
end

[Htot_ZF, ~] = freqz(conv(p, wZF.'), 1, 512);
plot(w/pi, 20*log10(abs(Htot_ZF)), 'b', 'LineWidth', 1.5); 
xlabel('Frecuencia normalizada (\times\pi rad/muestra)');
ylabel('dB');
legend('MMSE 25dB', 'MMSE 15dB', 'MMSE 5dB', 'ZF');
title('Respuestas en frecuencia - Canal global');
grid on;
