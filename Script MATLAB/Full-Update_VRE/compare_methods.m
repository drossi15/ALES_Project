%% COMPARE_METHODS.m
% Calcoliamo la PCA Statica (Modello fisso sui primi 200 campioni)
X_test = X_scaled(201:end, :); % Usiamo i dati già normalizzati
[P_stat, ~] = eig(corr(X_scaled(1:200, :))); 
[~, idx] = sort(diag(eig(corr(X_scaled(1:200, :)))), 'descend');
P_stat = P_stat(:, idx(1:n_pcs));

% Calcolo SPE statico per tutto il dataset di test
E_stat = X_test * (eye(m_vars) - P_stat * P_stat');
Q_stat_only = sum(E_stat.^2, 2);

% Plot di confronto
figure('Color', 'w', 'Name', 'Static vs Recursive PCA');
subplot(2,1,1);
plot(201:1000, Q_stat_only, 'r'); hold on;
yline(Q_lim_store(201), 'k--'); % Soglia fissa della PCA statica
% Aggiungiamo aree ombreggiate per i guasti
% Syntax: fill([x_start x_end x_end x_start], [y_min y_min y_max y_max], color)
yl = ylim;
patch([400 410 410 400], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([480 490 490 480], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([500 550 550 500], [yl(1) yl(1) yl(2) yl(2)], 'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([800 840 840 800], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title('STATIC PCA: Noticing False Alarms Due to Drift');
ylabel('SPE (Q)'); grid on;

subplot(2,1,2);
plot(201:1000, Q_store(201:1000), 'b'); hold on;
plot(201:1000, Q_lim_store(201:1000), 'r--');
% Aggiungiamo aree ombreggiate per i guasti
% Syntax: fill([x_start x_end x_end x_start], [y_min y_min y_max y_max], color)
yl = ylim;
patch([400 410 410 400], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([480 490 490 480], [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([500 550 550 500], [yl(1) yl(1) yl(2) yl(2)], 'y', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch([800 840 840 800], [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
title('RECURSIVE PCA: The model fits (No False Alarms)');
ylabel('SPE (Q)'); xlabel('Samples'); grid on;