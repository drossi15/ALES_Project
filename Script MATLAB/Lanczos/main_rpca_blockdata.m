%% MAIN_RPCA_ADAPTIVE_STRUCTURE.m
% Implementazione RPCA con Selezione Automatica del numero di PC (CPV)
% Prima occorre eseguire generate_data.m

if ~exist('X_raw', 'var')
    error('Esegui prima lo script generate_data.m!');
end


%% 1. Parametri dell'Algoritmo
mu = 0.95;              % Forgetting Factor (consiglio 0.95 per reattività)
alpha = 0.99;           % Confidenza per le soglie (99%)
cpv_threshold = 0.90;   % SOGLIA CPV: Vogliamo spiegare almeno il 90% della varianza
                        % Nota: Poiché il sistema ha 3 variabili latenti forti su 6, 
                        % il 90% dovrebbe selezionare n_pcs 


% Inizializzazione (Training Phase)
n_init = 200;
X_init = X_raw(1:n_init, :);
[n_obs, m_vars] = size(X_raw);

% Statistiche iniziali
b = mean(X_init)';      
sigma_vec = std(X_init)';
Sigma = diag(sigma_vec);
R = corr(X_init);

% Eigen-decomposition iniziale
[P_all, Lambda_all] = eig(R);
[lambda_sorted, idx] = sort(diag(Lambda_all), 'descend');
P_all = P_all(:, idx);





% SELEZIONE INIZIALE n_pcs (CPV)
total_variance = sum(lambda_sorted);
cum_variance = cumsum(lambda_sorted) / total_variance;
n_pcs = find(cum_variance >= cpv_threshold, 1, 'first');
if isempty(n_pcs), n_pcs = 1; end % Sicurezza

l_max = 5;

fprintf('Inizializzazione: Selezionate %d componenti principali (Spiegano %.2f%% var)\n', ...
    n_pcs, cum_variance(n_pcs)*100);

P = P_all(:, 1:n_pcs);



%% Loop Ricorsivo

B = 3;
buffer = [];

T2_store = zeros(n_obs, 1);
Q_store = zeros(n_obs, 1);
T2_lim_store = zeros(n_obs, 1);
Q_lim_store = zeros(n_obs, 1);
npcs_store = zeros(n_obs, 1); % Salviamo come cambia n_pcs nel tempo

fprintf('Avvio monitoraggio ricorsivo con struttura adattiva...\n');

for k = (n_init + 1) : n_obs
    
    % PASSO A: Acquisizione
    x_new_raw = X_raw(k, :)';

    % PASSO B: Detection (Usa il modello al passo k-1)
    x_new_scaled = (x_new_raw - b) ./ sigma_vec;

    % T2 Statistic
    scores = P' * x_new_scaled;
    Lambda_inv = diag(1./lambda_sorted(1:n_pcs));
    T2_stat = scores' * Lambda_inv * scores;
    
    % Q Statistic (SPE)
    x_hat = P * scores;
    residual = x_new_scaled - x_hat;
    Q_stat = residual' * residual;
    
    % Soglie Adattive
    T2_lim = chi2inv(alpha, n_pcs);
    
    % Soglia Q (Jackson-Mudholkar)
    theta1 = sum(lambda_sorted(n_pcs+1:end));
    theta2 = sum(lambda_sorted(n_pcs+1:end).^2);
    theta3 = sum(lambda_sorted(n_pcs+1:end).^3);
    h0 = 1 - (2*theta1*theta3)/(3*theta2^2);
    if h0 < 0, h0 = 0.001; end
    c_alpha = norminv(alpha);
    term1 = sqrt(2 * theta2 * h0^2) / theta1;
    term2 = 1 + (theta2 * h0 * (h0 - 1)) / (theta1^2);
    Q_lim = theta1 * (1 + term1 * c_alpha + term2)^(1/h0);
    
    % Salvataggio
    T2_store(k) = T2_stat;
    Q_store(k)  = Q_stat;
    T2_lim_store(k) = T2_lim;
    Q_lim_store(k)  = Q_lim;
    npcs_store(k) = n_pcs; % Salviamo il numero di PC usati
    
    % PASSO C: Logica di Aggiornamento
    is_fault = (T2_stat > T2_lim) || (Q_stat > Q_lim);
    
    if ~is_fault


        buffer = [buffer; x_new_raw'];
         
        
        if size(buffer,1) == B
            % 1. Aggiornamento Media e Varianza (Necessari per scalare i dati)
            b_old = b;
            b = mu * b_old + (1 - mu) * mean(buffer,1)';
            delta_b = b - b_old;
            
            sigma_sq_old = sigma_vec.^2;
            term_var = var(buffer,0,1)';
            sigma_sq_new = mu * (sigma_sq_old + delta_b.^2) + (1 - mu) * term_var;
            sigma_vec = sqrt(sigma_sq_new);
            
            % Aggiornamento R (Eq. 19)
            
            Sigma_new = diag(sigma_vec);
            Sigma_old = diag(sqrt(sigma_sq_old));
            inv_Sigma_new = diag(1 ./ sigma_vec);
            
            Term_Center = Sigma_old * R * Sigma_old + (delta_b * delta_b');
            Part1 = mu * (inv_Sigma_new * Term_Center * inv_Sigma_new);
            
            x_norm_update = (buffer - b') ./ sigma_vec'; 
            Part2 = (1 - mu) * (x_norm_update' * x_norm_update)/B;
            %disp(Part2)
    
            R = Part1 + Part2;
            R = (R + R') / 2; % Simmetria
            R = diag(1./sqrt(diag(R))) * R * diag(1./sqrt(diag(R))); %standardizzo per aver diagonale uguale a 1

            [P_all, Lambda_all] = eigs(R,l_max,'largestabs','IsSymmetric',true);
            [lambda_sorted, idx] = sort(diag(Lambda_all), 'descend');
            P_all = P_all(:, idx);
        
            % Calcolo CPV Ricorsivo
            total_var = trace(R); % somma degli elementi sulla diagonale
            cum_var_percentage = cumsum(lambda_sorted) / total_var;
           
            n_pcs_new = find(cum_var_percentage >= cpv_threshold, 1, 'first');
            n_pcs_new = min(n_pcs_new, l_max); % Non usare mai più di 4 PC su 6 sensori
            if isempty(n_pcs_new), n_pcs_new = l_max; end
            n_pcs = n_pcs_new;
            P = P_all(:, 1:n_pcs);

            buffer = [];    %Svuoto Buffer

        end
        
        
    else
        % Fault Detected: Modello congelato
    end
end




%% 3. Plotting Avanzato
figure('Name', 'Adaptive Structure RPCA', 'Color', 'w', 'Position', [100 100 1000 800]);

% Plot T2
subplot(3,1,1);
semilogy(T2_store, 'b'); hold on;
semilogy(T2_lim_store, 'r--');
% Aggiungiamo aree ombreggiate per i guasti
% Syntax: fill([x_start x_end x_end x_start], [y_min y_min y_max y_max], color)
yl = ylim;
patch([400 410 410 400], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([480 490 490 480], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([500 550 550 500], [yl(1) yl(1) yl(2) yl(2)], 'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([800 840 840 800], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title('T^2 Statistic (Variable Threshold)');
xlim([200 n_obs]); grid on;

% Plot Q
subplot(3,1,2);
semilogy(Q_store, 'b'); hold on;
semilogy(Q_lim_store, 'r--');
% Aggiungiamo aree ombreggiate per i guasti
% Syntax: fill([x_start x_end x_end x_start], [y_min y_min y_max y_max], color)
yl = ylim;
patch([400 410 410 400], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([480 490 490 480], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([500 550 550 500], [yl(1) yl(1) yl(2) yl(2)], 'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([800 840 840 800], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title('Q Statistic (SPE)');
xlim([200 n_obs]); grid on;

% Plot Numero di PC
subplot(3,1,3);
plot(npcs_store, 'k', 'LineWidth', 2);
ylim([0 m_vars]); 
yline(3, 'b:', 'Ground Truth (3)');
% Aggiungiamo aree ombreggiate per i guasti
% Syntax: fill([x_start x_end x_end x_start], [y_min y_min y_max y_max], color)
yl = ylim;
patch([400 410 410 400], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([480 490 490 480], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([500 550 550 500], [yl(1) yl(1) yl(2) yl(2)], 'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([800 840 840 800], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title(['Number of Principal Components (CPV Threshold = ' num2str(cpv_threshold*100) '%)']);
ylabel('n PCs'); xlabel('Samples');
xlim([200 n_obs]); grid on;
