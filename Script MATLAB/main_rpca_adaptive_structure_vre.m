

function [n_pcs_opt] = select_npcs_vre(R, P_all)

%--------------------------------------------------------------------------
% VRE-based selection of number of principal components
%
% INPUT:
%   R              : matrice di correlazione
%   P_all          : matrice degli Autovettori ordinati di R
%
% OUTPUT:
%   n_pcs_opt      : numero ottimo di pc (valore minimo di VRE)
%--------------------------------------------------------------------------

m = size(R,1);
VRE = zeros(m-1,1);
var_x = diag(R);

for l = 1:(m-1)

    P_l = P_all(:,1:l);    % PCA projection matrix
    C = P_l * P_l';        % C = P_l P_l^T
    u = zeros(m,1);        % varianza dell'errore di ricostruzione per ogni variabile

    for i = 1:m
        
        r_ii = R(i,i);
        r_i  = R(:,i);
        
        c_i  = C(:,i);
        c_ii = C(i,i);
        
        numeratore = r_ii - 2 * (c_i' * r_i) + (c_i' * R * c_i);
        denominatore = (1 - c_ii)^2;
        
        u(i) = numeratore / denominatore; % Eq(37)
    end
    VRE(l) = sum(u ./ var_x);  % Eq(38)
end

[~, n_pcs_opt] = min(VRE);  % il numero ottimo di pc è numero l che minimizza la VRE

end


%% MAIN_RPCA_ADAPTIVE_STRUCTURE.m

% Implementazione RPCA con Selezione Automatica del numero di PC con VRE
% Prima occorre eseguire generate_data.m

if ~exist('X_raw', 'var')
    error('Esegui prima lo script generate_data.m!');
end



%% 1. Parametri dell'Algoritmo
mu = 0.95;              % Forgetting Factor (consiglio 0.95 per reattività)
alpha = 0.99;           % Confidenza per le soglie (99%)

% Inizializzazione
n_init = 200;                    %%numero iniziale di campioni
X_init = X_raw(1:n_init, :);
[n_obs, m_vars] = size(X_raw);

%%

% Statistiche iniziali
b = mean(X_init)';            % Eq (3)   media delle 6 variabili
sigma_vec = std(X_init)';     %         std delle 6 variabili
Sigma = diag(sigma_vec);      % Eq (5)
R = corr(X_init);             % Eq(6)    matrice di correlazione

%%

% Scomposizione Autovalori e Autovettori  Eq(1) e Eq(2)

[P_all, Lambda_all] = eig(R);    %P_all matrice le cui colonne autovettori di R, lambda matrice diagonale degli autovalori   
[lambda_sorted, idx] = sort(diag(Lambda_all), 'descend');    %Estraggo gli autovalori e li riordino dal più grande al più piccolo
P_all = P_all(:, idx);   %Riordino le colonne di P_all usando idx ( j-esima colonna di P_all corrisponde a lambda_sorted(j))


%% --- SELEZIONE INIZIALE n_pcs (VRE) ---   Eq(33)


n_pcs= select_npcs_vre(R, P_all); %Selezioni il numero di componenti con VRE
P = P_all(:, 1:n_pcs);   % Eq(1)  %Selezione le prime n_pcs colonne di P, gli autovettori corrispondenti alle prime n_pcs componenti principali




%% 2. Loop Ricorsivo (Adaptive Monitoring & Structure)
T2_store = zeros(n_obs, 1);
Q_store = zeros(n_obs, 1);
T2_lim_store = zeros(n_obs, 1);
Q_lim_store = zeros(n_obs, 1);
npcs_store = zeros(n_obs, 1); % Salviamo come cambia n_pcs nel tempo

fprintf('Avvio monitoraggio ricorsivo con struttura adattiva...\n');

for k = (n_init + 1) : n_obs
    
    % Acquisizione nuova osservazione
    x_new_raw = X_raw(k, :)';
    
    % Eseguo lo Scaling del dato
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
    theta1 = sum(lambda_sorted(n_pcs+1:end));     % Eq(35)
    theta2 = sum(lambda_sorted(n_pcs+1:end).^2);  % Eq(36)
    theta3 = sum(lambda_sorted(n_pcs+1:end).^3);  % Eq(42)
    h0 = 1 - (2*theta1*theta3)/(3*theta2^2);      % Eq(41)
    if h0 < 0, h0 = 0.001; end
    c_alpha = norminv(alpha);
    term1 = sqrt(2 * theta2 * h0^2) / theta1;
    term2 = 1 + (theta2 * h0 * (h0 - 1)) / (theta1^2);
    Q_lim = theta1 * (1 + term1 * c_alpha + term2)^(1/h0);  % Eq(40)
    
    % Salvataggio
    T2_store(k) = T2_stat;
    Q_store(k)  = Q_stat;
    T2_lim_store(k) = T2_lim;
    Q_lim_store(k)  = Q_lim;
    npcs_store(k) = n_pcs;  %Salvo il numero di pc usato
    
    % Logica di Aggiornamento %step 3 section 7.2
    is_fault = (T2_stat > T2_lim) || (Q_stat > Q_lim);  
    
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
       
        n_pcs= select_npcs_vre(R, P_all);
        P = P_all(:, 1:n_pcs);
        
    else
        % Fault Detected: Stop Update
        % n_pcs rimane quello dell'ultimo passo sano
        % Fault Detected: Log the fault occurrence
        fprintf('Fault detected at observation %d\n', k);
    end
end



%% 3. Plotting Avanzato
figure('Name', 'Adaptive Structure RPCA', 'Color', 'w', 'Position', [100 100 1000 800]);

% Plot T2
subplot(3,1,1);
semilogy(T2_store, 'b'); hold on;
semilogy(T2_lim_store, 'r--');
xline(800, 'k-');
title('T^2 Statistic (Variable Threshold)');
xlim([200 n_obs]); grid on;

% Plot Q
subplot(3,1,2);
semilogy(Q_store, 'b'); hold on;
semilogy(Q_lim_store, 'r--');
xline(800, 'k-');
title('Q Statistic (SPE)');
xlim([200 n_obs]); grid on;

% Plot Numero di PC
subplot(3,1,3);
plot(npcs_store, 'k', 'LineWidth', 2);
ylim([0 m_vars]); 
yline(3, 'b:', 'Ground Truth (3)');
xline(800, 'r-', 'Fault Start');
title(['Number of Principal Components (CPV Threshold']);
ylabel('n PCs'); xlabel('Samples');
xlim([200 n_obs]); grid on;