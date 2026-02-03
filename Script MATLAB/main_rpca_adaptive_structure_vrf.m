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
                        % il 90% dovrebbe selezionare n_pcs = 3.

% -- Inizializzazione (Training Phase) --
n_init = 200;
X_init = X_raw(1:n_init, :);
[n_obs, m_vars] = size(X_raw);

%%

% Statistiche iniziali
b = mean(X_init)';         %Eq (3)   
sigma_vec = std(X_init)';   
Sigma = diag(sigma_vec);   %Eq (5)
R = corr(X_init);         %Eq(6) Correlation Matrix

% Eigen-decomposition iniziale    Eq(1) e Eq(2)
[P_all, Lambda_all] = eig(R);
[lambda_sorted, idx] = sort(diag(Lambda_all), 'descend');
P_all = P_all(:, idx);



% --- SELEZIONE INIZIALE n_pcs (CPV) ---   Eq(33)
total_variance = sum(lambda_sorted);
cum_variance = cumsum(lambda_sorted) / total_variance;
n_pcs = find(cum_variance >= cpv_threshold, 1, 'first');
if isempty(n_pcs), n_pcs = 1; end % Sicurezza

fprintf('Inizializzazione: Selezionate %d componenti principali (Spiegano %.2f%% var)\n', ...
    n_pcs, cum_variance(n_pcs)*100);

P = P_all(:, 1:n_pcs);   % Eq(1)