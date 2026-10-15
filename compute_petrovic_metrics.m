function [Q_ABF, L_ABF, N_ABF, MAD_val, SD_val] = compute_petrovic_metrics(imgA, imgB, imgF)
% 1. Costanti sigmoidali standard di Petrovic & Xydeas (2005)
L = 1.5;
Gamma_g = 0.9994;  kappa_g = 15;  sigma_g = 0.5;
Gamma_a = 0.9879;  kappa_a = 22;  sigma_a = 0.8;

% 2. Operatori di gradiente di Sobel 3x3
h_x = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
h_y = h_x';

% 3. Calcolo delle derivate spaziali con i filtri di Sobel
sx_A = imfilter(imgA, h_x, 'replicate');  sy_A = imfilter(imgA, h_y, 'replicate');
sx_B = imfilter(imgB, h_x, 'replicate');  sy_B = imfilter(imgB, h_y, 'replicate');
sx_F = imfilter(imgF, h_x, 'replicate');  sy_F = imfilter(imgF, h_y, 'replicate');

% 4. Modulo del gradiente (forza dei bordi g)
g_A = sqrt(sx_A.^2 + sy_A.^2);
g_B = sqrt(sx_B.^2 + sy_B.^2);
g_F = sqrt(sx_F.^2 + sy_F.^2);

% 5. Orientazione angolare dei bordi (alpha)
a_A = atan2(sy_A, sx_A);
a_B = atan2(sy_B, sx_B);
a_F = atan2(sy_F, sx_F);

%atan2 anziché atan (la funzione per l'arcotangente) è numericamente
%robusta e calcola l'angolo corretto su
%tutti e 4 i quadranti ed evita la divisione per zero quando s_x=0.

% 6. Calcolo della conservazione dei bordi per l'immagine A (Q_AF)
% A. Rapporto relativo di forza del gradiente G_AF
max_g_AF = max(g_A, g_F);
min_g_AF = min(g_A, g_F);
G_AF = zeros(size(g_A));
val_AF = max_g_AF > 0;
G_AF(val_AF) = min_g_AF(val_AF) ./ max_g_AF(val_AF);

% B. Differenza angolare relativa Delta_A_AF (con simmetria a pi greco)
diff_a_AF = abs(a_A - a_F);
diff_a_AF(diff_a_AF > pi/2) = pi - diff_a_AF(diff_a_AF > pi/2);
Delta_A_AF = 1 - (diff_a_AF / (pi / 2));

% C. Sigmoidi percettive e punteggio totale Q_AF
Qg_AF = Gamma_g ./ (1 + exp(-kappa_g * (G_AF - sigma_g)));
Qa_AF = Gamma_a ./ (1 + exp(-kappa_a * (Delta_A_AF - sigma_a)));
Q_AF  = Qg_AF .* Qa_AF;

% 7. Calcolo della conservazione dei bordi per l'immagine B (Q_BF)
max_g_BF = max(g_B, g_F);
min_g_BF = min(g_B, g_F);
G_BF = zeros(size(g_B));
val_BF = max_g_BF > 0;
G_BF(val_BF) = min_g_BF(val_BF) ./ max_g_BF(val_BF);

diff_a_BF = abs(a_B - a_F);
diff_a_BF(diff_a_BF > pi/2) = pi - diff_a_BF(diff_a_BF > pi/2);
Delta_A_BF = 1 - (diff_a_BF / (pi / 2));

Qg_BF = Gamma_g ./ (1 + exp(-kappa_g * (G_BF - sigma_g)));
Qa_BF = Gamma_a ./ (1 + exp(-kappa_a * (Delta_A_BF - sigma_a)));
Q_BF  = Qg_BF .* Qa_BF;

% 8. Pesi percettivi di salienza visiva w_A e w_B (potenza L = 1.5)
w_A = g_A .^ L;
w_B = g_B .^ L;

% 9. Denominatore comune (somma totale dei pesi percettivi)
denom = sum(w_A(:) + w_B(:));
if denom == 0, denom = 1e-12; end

% 10. Calcolo di Q^{AB/F} (Fusion Information Score, Eq. 11)
Q_ABF = sum(Q_AF(:) .* w_A(:) + Q_BF(:) .* w_B(:)) / denom;

% 11. Calcolo di L^{AB/F} (Fusion Loss, Eq. 12-13)
% Flag indicatore r: vale 1 dove il gradiente fuso si è attenuato rispetto a uno degli input
r = double((g_F < g_A) | (g_F < g_B));
loss_local = (1 - Q_AF) .* w_A + (1 - Q_BF) .* w_B;
L_ABF = sum(r(:) .* loss_local(:)) / denom;

% 12. Calcolo di N^{AB/F} (Fusion Artifacts / Noise, Eq. 14-15)
% Indicatore di artefatto: attivo dove il gradiente fuso è strettamente MAGGIORE di ENTRAMBI gli input
art_mask = (g_F > g_A) & (g_F > g_B);
N_local = zeros(size(g_F));
N_local(art_mask) = 2 - Q_AF(art_mask) - Q_BF(art_mask);
N_ABF = sum(N_local(:) .* (w_A(:) + w_B(:))) / denom;

% 13. Metriche Statistiche: MAD (Eq. 16) e SD (Eq. 17)
MAD_val = mean(abs(imgA(:) - imgF(:)) + abs(imgB(:) - imgF(:))) / 2;
SD_val  = std(imgF(:));
end
