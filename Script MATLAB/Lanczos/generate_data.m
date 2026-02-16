%% 1. Configurazione Iniziale
clear; clc; close all;
rng(158);

% Parametri generali
n_samples = 1000;      % Numero totale di campioni
n_vars = 6;            % Numero di sensori (x1...x6)
n_latents = 3;         % Numero di variabili latenti (t1...t3)

% Matrice di Mixing A (Dalla Slide 6 del Project Guidelines)
A_nominal = [
   -0.2310, -0.0816, -0.2662;
   -0.3241,  0.7055, -0.2158;
   -0.2170, -0.3056, -0.5207;
   -0.4089, -0.3442, -0.4501;
   -0.6408,  0.3102,  0.2372;
   -0.4655, -0.4330,  0.5938
];

% Deviazioni standard delle variabili latenti t (Dalla Slide 6)
std_t = [1.0; 0.8; 0.6]; 
% Deviazione standard del rumore (Dalla Slide 6)
std_noise = 0.2;

%% 2. Generazione Dati Time-Varying
% Inizializziamo la matrice dei dati X (Righe=Campioni, Colonne=Variabili)
X_raw = zeros(n_samples, n_vars);
GroundTruth_Fault = zeros(n_samples, 1); % Per verifica finale (1=Guasto, 0=Sano)

% Parametri per la simulazione "Time-Varying" 
% Simuliamo un invecchiamento lento del processo.
% Modifichiamo lentamente l'elemento (1,1) della matrice A
drift_slope_A = 0.0002; 
% Simuliamo una deriva del sensore 2 (spostamento della media)
sensor_drift_slope = 0.002; 


fprintf('Generazione dati in corso...\n');

for k = 1:n_samples
    
    %A. Variabili Latenti (t)
    % Generiamo t1, t2, t3 con le deviazioni standard specificate
    t = randn(n_latents, 1) .* std_t;
    
    % B. Rendere il sistema Time-Varying (Richiesta Progetto)
    % 1. Modifica della struttura di correlazione (Matrice A cambia piano piano)
    A_current = A_nominal;
    A_current(1,1) = A_nominal(1,1) + (k * drift_slope_A); 
    
    % C. Calcolo delle uscite x (senza rumore)
    x_k = A_current * t;
    
    % 2. Aggiunta deriva sulla media (Sensor Drift su x2)
    % Questo simula lo sporcamento del sensore citato nel paper
    x_k(2) = x_k(2) + (k * sensor_drift_slope);
    
    % D. Aggiunta Rumore
    noise = randn(n_vars, 1) * std_noise;
    x_k = x_k + noise;
    
    % E. Inserimento Guasto (Fault Injection)
    
    % CASO 1: Intermittent Faults (Spikes)
    % Simulazione di disturbi impulsivi su Sensore 5 (es. interferenze elettriche)
    if k >= 400 && k<=410 || (k >= 480 && k <= 490)
        x_k(5) = x_k(5) + 5; % Magnitude elevata ma breve durata
        GroundTruth_Fault(k) = 1;
    end
    
    % CASO 2: Incipient Fault (Drift temporaneo / "Dome Shape")
    % Simulazione di un'anomalia che cresce e poi rientra (es. surriscaldamento temporaneo)
    if k >= 500 && k <= 550
        % Crea una curva che va da 0 a 6.0 e torna a 0 tra k=500 e k=550
        duration = 550 - 500;
        t_local = k - 500;
        fault_magnitude = 6.0 * sin(pi * t_local / duration); 
        
        x_k(4) = x_k(4) + fault_magnitude; 
        GroundTruth_Fault(k) = 1;
    end

    % CASO 3: Permanent Fault (Step Change)
    % Simulazione di rottura o bias permanente sul Sensore 1
    if k >= 800 && k<=840
        x_k(1) = x_k(1) + 4.0;
        GroundTruth_Fault(k) = 1;
    end
    
    % Salvataggio
    X_raw(k, :) = x_k';
end

%% 3. Pre-processing Iniziale (Scaling)

% Usiamo i primi 200 campioni come "Training Set" iniziale per scalare
n_train = 200;
X_train = X_raw(1:n_train, :);

% Calcoliamo media e std sul training set (sano e iniziale)
mean_vec = mean(X_train);
std_vec = std(X_train);

% Normalizziamo TUTTO il dataset usando le statistiche iniziali
X_scaled = (X_raw - mean_vec) ./ std_vec;

%% 4. Visualizzazione
figure('Name', 'Simulated Data Time-Varying', 'Color', 'w');
subplot(2,1,1);
plot(X_scaled);
title('Simulated Process Data (Normalized)');
ylabel('Normalized Amplitude');
xlabel('Samples (k)');
grid on;

subplot(2,1,2);
% Zoom sul sensore 2 per vedere il drift
plot(X_scaled(:,2), 'b'); hold on;
yline(0, 'k--');
title('Zoom Sensor 2 (Note the slow upward Drift)');
xlabel('Samples (k)');
grid on;

fprintf('Dati generati. Matrice X_scaled pronta (%dx%d).\n', size(X_scaled));