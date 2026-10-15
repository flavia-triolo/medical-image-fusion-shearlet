function W = compute_gfs_weights(base_layer, m)  %prende in ingresso il base layer 256x256 e la dimensione m=5
% 1. Dimensioni dell'immagine e raggio della finestra locale
[H, W_dim] = size(base_layer); % restituisce il numero di righe e di colonne (H=256, W_dim=256)
r = floor(m / 2); 
% se m=5, la funzionione (approssimazione per difetto) restituisce 2 (intorno -2, +2 del pixel centrale)

% 2. Padding simmetrico ai bordi per evitare errori sui margini
%img_padded aggiunge r=2 righe sopra, 2 righe sotto, 2 colonne a sinistra e
%2 colonne a destra. La matrice allargata diventa d dimensione
%(256+4)x(256+4). L'opzione 'symmetric' riempie la cornice specchiando i
%pixel del bordo così da non creare salti artificiali o cornici nere che
%altererebbero il calcolo della covarianza.
img_padded = padarray(base_layer, [r, r], 'symmetric');  
W = zeros(H, W_dim);

% 3. Scansione di tutti i pixel dell'immagine
% due cicli annidati che scorrono tutte le 256 righe e tutte le 256 colonne
% dell'immagine.
    for i = 1:H
        for j = 1:W_dim
            % Estrazione della finestra locale 5x5 centrata in (i, j)
            Z = img_padded(i : i + 2*r, j : j + 2*r);
            % se r=2, l'intervallo i:i+4 seleziona esattamente 5 righe e
            % j:j+4 seleziona 5 colonne
            % della matrice con padding. Z è la matrice 5x5 centrata sul
            % pixel (i,j).
            
            % Forza dei bordi orizzontali alpha_H (Eq. 4-6)
            % (Le righe sono le osservazioni, le colonne sono le variabili)
            C_H = cov(Z);
            % alpha_H = sum(eig(C_H));
            % in algebra la somma degli autovalori di una matrice quadrata
            % coincide con la sua traccia e riduce il costo computazionale.
            alpha_H = trace(C_H);
            % Forza dei bordi verticali alpha_V (Eq. 7)
            % (Le colonne sono le osservazioni, le righe sono le variabili -> trasposta Z')
            C_V = cov(Z');
            alpha_V = trace(C_V);
            
            % Peso statistico complessivo per il pixel (i, j) (Eq. 8)
            W(i, j) = alpha_H + alpha_V;
        end
    end
end