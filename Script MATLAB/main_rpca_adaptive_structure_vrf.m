%% MAIN_RPCA_ADAPTIVE_STRUCTURE.m
% Implementazione RPCA con Selezione Automatica del numero di PC con VRF
% Prima occorre eseguire generate_data.m

if ~exist('X_raw', 'var')
    error('Esegui prima lo script generate_data.m!');
end

%% 1. Parametri dell'Algoritmo
mu = 0.95;              % Forgetting Factor (consiglio 0.95 per reattività)
alpha = 0.99;           % Confidenza per le soglie (99%)
cpv_threshold = 0.90;   % SOGLIA CPV: Vogliamo spiegare almeno il 90% della varianza
                        % Nota: Poiché il sistema ha 3 variabili latenti forti su 6, 
                        % il 90% dovrebbe selezionare n_pcs = 3.

% Inizializzazione
n_init = 200;                    %%numero iniziale di campioni
X_init = X_raw(1:n_init, :);
[n_obs, m_vars] = size(X_raw);

%%

% Statistiche iniziali
b = mean(X_init)';            %Eq (3)   media delle 6 variabili
sigma_vec = std(X_init)';     %         std delle 6 variabili
Sigma = diag(sigma_vec);      %Eq (5)
R = corr(X_init);             %Eq(6)    matrice di correlazione

%%

% Scomposizione Autovalori e Autovettori  Eq(1) e Eq(2)

[P_all, Lambda_all] = eig(R);    %P_all matrice le cui colonne autovettori di R, lambda matrice diagonale degli autovalori   
[lambda_sorted, idx] = sort(diag(Lambda_all), 'descend');    %Estraggo gli autovalori e li riordino dal più grande al più piccolo
P_all = P_all(:, idx);   %Riordino le colonne di P_all usando idx ( j-esima colonna di P_all corrisponde a lambda_sorted(j))


%% --- SELEZIONE INIZIALE n_pcs (CPV) ---   Eq(33)

total_variance = sum(lambda_sorted);     % varianza totale come somma degli autovalori (teoricamente 6)
cum_variance = cumsum(lambda_sorted) / total_variance;  %frazioni cumulative (vettore 6*1)
n_pcs = find(cum_variance >= cpv_threshold, 1, 'first');   %trova il numero minimo numero di componenti per superare la soglia cpv
if isempty(n_pcs), n_pcs = 1; end % Sicurezza

fprintf('Inizializzazione: Selezionate %d componenti principali (Spiegano %.2f%% var)\n', ...
    n_pcs, cum_variance(n_pcs)*100);

P = P_all(:, 1:n_pcs);   % Eq(1)  %Selezione le prime n_pcs colonne di P, gli autovettori corrispondenti alle prime n_pcs componenti principali


%%  Selezione iniziale con VRF

u = zeros(m_vars);

 for i = 1: m_vars:
    u(i)= Lambda_all(i)- 
 end




%% 2. Loop Ricorsivo (Adaptive Monitoring & Structure)
T2_store = zeros(n_obs, 1);
Q_store = zeros(n_obs, 1);
T2_lim_store = zeros(n_obs, 1);
Q_lim_store = zeros(n_obs, 1);
npcs_store = zeros(n_obs, 1); % Salviamo come cambia n_pcs nel tempo

fprintf('Avvio monitoraggio ricorsivo con struttura adattiva...\n');

for k = (n_init + 1) : n_obs
    
    %PASSO A: Acquisizione nuova osservazione
    x_new_raw = X_raw(k, :)';
    
    % --- PASSO B: Detection (Usa il modello al passo k-1) ---
    x_new_scaled = (x_new_raw - b) ./ sigma_vec;
    
    % T2 Statistic     Eq (43)
    scores = P' * x_new_scaled;
    Lambda_inv = diag(1./lambda_sorted(1:n_pcs));
    T2_stat = scores' * Lambda_inv * scores;
    
    % Q Statistic (SPE)   Eq (39)
    x_hat = P * scores;
    residual = x_new_scaled - x_hat;
    Q_stat = residual' * residual;
    
    % Soglie Adattive, section 7.1
    T2_lim = chi2inv(alpha, n_pcs);   
    
    % Soglia Q (Jackson-Mudholkar)
    theta1 = sum(lambda_sorted(n_pcs+1:end));    %Eq(35)
    theta2 = sum(lambda_sorted(n_pcs+1:end).^2);  %Eq(36)
    theta3 = sum(lambda_sorted(n_pcs+1:end).^3);  %Eq(42)
    h0 = 1 - (2*theta1*theta3)/(3*theta2^2);    %Eq(41)
    if h0 < 0, h0 = 0.001; end
    c_alpha = norminv(alpha);
    term1 = sqrt(2 * theta2 * h0^2) / theta1;
    term2 = 1 + (theta2 * h0 * (h0 - 1)) / (theta1^2);
    Q_lim = theta1 * (1 + term1 * c_alpha + term2)^(1/h0); %Eq(40)
    
    % Salvataggio
    T2_store(k) = T2_stat;
    Q_store(k)  = Q_stat;
    T2_lim_store(k) = T2_lim;
    Q_lim_store(k)  = Q_lim;
    npcs_store(k) = n_pcs; % Salviamo il numero di PC usati
    
    % --- PASSO C: Logica di Aggiornamento ---
    is_fault = (T2_stat > T2_lim) || (Q_stat > Q_lim);  %step 3 section 7.2
    
    if ~is_fault
        % Aggiornamento Media e Varianza
        b_old = b;
        b = mu * b_old + (1 - mu) * x_new_raw; %Eq (15)
        delta_b = b - b_old;
        
        sigma_sq_old = sigma_vec.^2;
        term_var = (x_new_raw - b).^2;
        sigma_sq_new = mu * (sigma_sq_old + delta_b.^2) + (1 - mu) * term_var; %Eq(16)
        sigma_vec = sqrt(sigma_sq_new);
        
        % Aggiornamento R (Eq. 19)
        Sigma_new = diag(sigma_vec);
        Sigma_old = diag(sqrt(sigma_sq_old));
        inv_Sigma_new = diag(1 ./ sigma_vec);
        
        Term_Center = Sigma_old * R * Sigma_old + (delta_b * delta_b');
        Part1 = mu * (inv_Sigma_new * Term_Center * inv_Sigma_new);
        
        x_norm_update = (x_new_raw - b) ./ sigma_vec; 
        Part2 = (1 - mu) * (x_norm_update * x_norm_update');
        
        R = Part1 + Part2;
        R = (R + R') / 2; % Simmetria
        
        % --- EIGEN DECOMPOSITION & SELEZIONE RICORSIVA n_pcs ---
        [P_all, Lambda_all] = eig(R);
        [lambda_sorted, idx] = sort(diag(Lambda_all), 'descend');
        P_all = P_all(:, idx);
        
        % Calcolo CPV Ricorsivo
        total_var = sum(lambda_sorted);
        cum_var_percentage = cumsum(lambda_sorted) / total_var;
        
        % Aggiorniamo n_pcs per il prossimo ciclo
        n_pcs_new = find(cum_var_percentage >= cpv_threshold, 1, 'first');
        if isempty(n_pcs_new), n_pcs_new = 1; end
        
        % Aggiorniamo le variabili globali del modello
        n_pcs = n_pcs_new;
        P = P_all(:, 1:n_pcs);
        
    else
        % Fault Detected: Stop Update
        % n_pcs rimane quello dell'ultimo passo sano
    end
end