function [P_new, Lambda_new] = rank_one_update(P_old, Lambda_old, z, mu)
    % RANK_ONE_UPDATE Implementazione dell'aggiornamento ricorsivo (Li et al.)
    % Input:
    %   P_old: Matrice autovettori precedente (m x m)
    %   Lambda_old: Matrice diagonale autovalori precedente (m x m)
    %   z: Nuovo campione normalizzato (colonna m x 1)
    %   mu: Forgetting factor
    
    % 1. Proiettiamo il nuovo dato nello spazio degli autovettori attuali
    % Questo ci dice "quanto" il nuovo dato sposta ogni componente.
    v = P_old' * z; 
    
    % 2. Definiamo il peso del nuovo dato
    gamma = 1 - mu;
    
    % 3. Costruiamo la matrice "Diagonal plus Rank-One"
    % Invece di fare eig(R) che è una matrice densa, facciamo eig su M.
    % M è diagonale (mu * Lambda_old) + una piccola correzione (gamma * v * v').
    % Calcolare gli autovalori di una matrice quasi diagonale è molto più veloce.
    M = mu * Lambda_old + gamma * (v * v');
    
    % 4. Calcoliamo la decomposizione di M
    [P_rot, Lambda_new] = eig(M);
    
    % 5. Ruotiamo i vecchi autovettori per ottenere quelli nuovi
    P_new = P_old * P_rot;
    
    % Nota: MATLAB non garantisce l'ordine decrescente, quindi 
    % l'ordinamento lo gestiremo nel Main per sicurezza.

end