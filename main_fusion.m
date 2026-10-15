%Questo codice crea in output una cartella in cui vengono salvate le 15
%immagini fuse. Inoltre si calcolano le metriche sull'intero dataset fuso.

%In generale all'inizio di ogni script MATLAB si inseriscono i seguenti tre
%comandi per assicurarsi che la memoria sia pulita e non ci siano residui
%di esecuzioni precedenti:

clc;         %sta per Clear Command Window: cancella tutto il testo visibile nella finestra dei comandi in basso, rendendola pulita
clear;       %cancella tutte le variabili salvate nella memoria di lavoro (workspace). Evita bug (ad esempio riutilizzare per sbaglio una variabile definita in un calcolo precedente)
close all;   %chiude tutte le eventuali finestre di immagini o grafici rimaste aperte

%NOTA: nella sintassi MATLAB il ";" serve a non stampare il risultato a
%schermi nella Command Window.
%Adesso si deve dire a MATLAB dove si trovano le cartelle con le immagini e
%dove si vogliono salvare i risultati.

ct_dir = 'Medical-image-fusion-dataset-master/test/CT';
mri_dir = 'Medical-image-fusion-dataset-master/test/MRI';
out_dir = 'output_fused';

%NOTA: in MATLAB non serve dichiarare il tipo di variabile (come int,
%string o float in C). Si sta memorizzando inoltre il percorso relativo
%della cartella delle immagini che parte dalla posizione in cui si trova lo
%script.


if ~exist(out_dir, 'dir')
    mkdir(out_dir);       % Crea la cartella di output se non esiste.
end


%exist(nome, 'tipo') è una funzione MATLAB che controlla se un file o una
%cartella esiste. In MATLAB il simbolo ~ è l'operatore logico NOT
%Adesso si deve far leggere il contenuto della cartella e contare quanti
%file ci sono, memorizzandone l'elenco.

ct_files = dir(fullfile(ct_dir, '*.jpg'));  %* fa in modo che vengano guardati tutti i file che terminano con l'estensione .jpg
mri_files = dir(fullfile(mri_dir, '*.jpg'));
num_images = length(ct_files);  %calcola il numero di elementi che viene salvato dentro la variabile num_images
% fullfile è una funzione che unisce cartelle e nomi file.
%dir è il comando che interroga il sistema operativo e trova tutti i file
%corrispondenti che salva il risultato un una struttura struct array, che
%contiene per ogni file il nome, la dimensione, la data di modifica, etc..

%Definisco i parametri dell'algoritmo:

nScales = 4;        %imposta il numero di scale di decomposizione a J=4
shearLevels = [1, 1, 2, 2]; %è un vettore riga. Alle prime due scale si vuole shear level 1 (2^{1+2}=8 direzioni), e alle ultime due scale shear level 2

r_gf = 30;     %raggio per il guided filter
eps_gf = 0.1;  %parametri di regolarizzazione epsilon per il guided filter
m_gfs = 5;     %dimensione della finestra per la statistica GFS


saliency_on_base = true; % false: calcola salienza su imgA/imgB (Input); true: calcola su baseA/baseB (Base layer)
base_fusion_mode = 2; % 1: Formulazione a Pesi Congiunti; 2: Formulazione a Due Stadi

% Costruzione del sistema di filtri Shearlet (griglia 256x256)
shearletSystem = SLgetShearletSystem2D(0, 256, 256, nScales, shearLevels);

% Matrice per salvare le 5 metriche per tutte le 15 immagini [Q, L, N, MAD, SD]
metrics_results = zeros(num_images, 5);
for k = 1:num_images
    fprintf('Elaborazione coppia %d di %d: %s\n', k, num_images, ct_files(k).name);

    imgA = imread(fullfile(ct_dir, ct_files(k).name));
    imgB = imread(fullfile(mri_dir, mri_files(k).name));

    if size(imgA, 3) == 3, imgA = rgb2gray(imgA); end   %size(imgA, 3) controlla la terza dimensione della matrice: se è uguale a 3, significa che l'immagine ha 3 canali di colore, allora la funzione rgb2gray la trasforma in una matrice bidimensionale a intensità di grigi.
    if size(imgB, 3) == 3, imgB = rgb2gray(imgB); end
    
    %Converto i valori da interi in reali (double), mappando l'intervallo
    %[0, 255] nell'intervallo reale [0, 1]. (I_double=I_unit8/255)
    imgA = im2double(imgA);
    imgB = im2double(imgB);

    % Calcolo dei coefficienti Shearlet per la coppia k-esima (eseguito 15 volte!)
    coeffsA = SLsheardec2D(imgA, shearletSystem);
    coeffsB = SLsheardec2D(imgB, shearletSystem);

    % 1. Estrazione del Base Layer (bassa frequenza, ULTIMA fetta)
    % NOTA IMPORTANTE: in questa configurazione di ShearLab il base layer
    % (passa-basso) è memorizzato nell'ULTIMA fetta del tensore, non nella
    % prima come si potrebbe pensare per analogia con la wavelet.
    % Verificato empiricamente: la fetta a energia massima è la 49/49
    % (con nScales=4, shearLevels=[1,1,2,2]), non la 1/49.
    idx_base = size(coeffsA, 3); % calcolato dinamicamente, non hardcodato a 49
    baseA = coeffsA(:, :, idx_base);
    baseB = coeffsB(:, :, idx_base);

    % 2. Inizializzazione del tensore per i coefficienti fusi
    coeffsF = zeros(size(coeffsA));
    % size(coeffsA)
    % restituisce le dimensioni del tensore (cioè 256x256x49).
    % zeros(...) alloca in memoria una nuova matrice 3D piena di zeri con
    % la stessa identica dimensione, pronta ad accogliere:
    % nell'ultima fetta (idx_base) il base layer unificato fuso (L_unified)
    % nelle fette da 1 a idx_base-1 le sotto-bande di dettaglio fuse con choose-max
    % 3. Kernel binomiale 5x5 separabile (Achanta et al., 2009)
    h_binom = [1, 4, 6, 4, 1]' * [1, 4, 6, 4, 1] / 256;

    % Scelta della sorgente per la salienza (Input originale vs Base layer)

    if saliency_on_base
        sal_inA = baseA;
        sal_inB = baseB;
    else
        sal_inA = imgA;
        sal_inB = imgB;
    end

    % 4. Calcolo delle mappe di salienza S_A e S_B (Eq. 1)
    S_A = abs(mean(sal_inA(:)) - imfilter(sal_inA, h_binom, 'replicate', 'same'));
    S_B = abs(mean(sal_inB(:)) - imfilter(sal_inB, h_binom, 'replicate', 'same'));

    %mean(...) calcola la media globale scalare I_mu
    %imfilter esegue la convoluzione 2D spaziale dell'immagine con il
    %filtro H: 'replicate' gestisce i bordi dell'immagine replicando i
    %pixel periferici (evitando che si formino cornici scure artificiali ai
    %margini), 'same' garantisce che l'immagine filtrata I_f abbia
    %esattamente le stesse dimensioni dell'originale, abs calcola il valore
    %assoluto elemento per elemento (la norma L_2 per scalari).


    % 5. Generazione delle mappe binarie (Eq. 2)
    BW_A = double(S_A >= S_B);
    BW_B = 1 - BW_A;
    % S_A >= S_B esegue un confronto elemento per elemento tra tutti i
    % pixel delle due matrici e restituisce una matrice booleana di tipo
    % logical contenente true dove la salienza di A è maggiore o uguale a
    % quella di B, e false altrove. double(...) trasforma i valori logici
    % true/false nei numeri reali 1 e 0, necessari per le operazioni
    % matematiche successive.

    % 6. Regolarizzazione spaziale tramite Guided Filter (Eq. 3)
    W_A = guided_filter(sal_inA, BW_A, r_gf, eps_gf);
    W_B = guided_filter(sal_inB, BW_B, r_gf, eps_gf);
    % BW_A e BW_B sono gli input da filtrare poichè presentano transizioni brusche 
    % e rumore locale. Il guided filter usa la
    % struttura dell'immagine guida per lisciare la mappa; sal_inA e
    % sal_inB sono le guide.

    % Alternativa nativa con Image Processing Toolbox:
    % q = imguidedfilter(p, I, 'NeighborhoodSize', [2*r+1, 2*r+1], 'DegreeOfSmoothing', eps_val);

    % 7. Calcolo dei pesi statistici GFS (matrici di covarianza ed autovalori, Eq. 4-8)
    W_stat_A = compute_gfs_weights(baseA, m_gfs);
    W_stat_B = compute_gfs_weights(baseB, m_gfs);

    % 8. Fusione del Base Layer (Eq. 9)
    delta = 1e-8; % Costante di stabilizzazione per evitare divisioni per zero

    if base_fusion_mode == 1
        % Formulazione 1: Pesi congiunti (Prodotto Salienza-GF e Statistica GFS)
        W_tot_A = W_A .* W_stat_A;
        W_tot_B = W_B .* W_stat_B;
        base_fused = (baseA .* W_tot_A + baseB .* W_tot_B) ./ (W_tot_A + W_tot_B + delta);
    else
        % Formulazione 2: Fusione a due stadi (Pseudocodice pag. 7)
        F_gfs = (baseA .* W_stat_A + baseB .* W_stat_B) ./ (W_stat_A + W_stat_B + delta);
        F_sal = (baseA .* W_A + baseB .* W_B) ./ (W_A + W_B + delta);
        base_fused = 0.5 * F_gfs + 0.5 * F_sal;
    end

    % Salvataggio del Base Layer fuso nell'ultima fetta di coeffsF
    coeffsF(:, :, idx_base) = base_fused;

    % i coefficienti di dettaglio di una shearlet possono essere sia
    % positivi che
    % negativi, in
    % quanto raprresentano la singolarità direzionale del fronte d'onda.
    % L'intensità o la forza
    % di un bordo o di
    % dettaglio anatomico è misurata dal suo valore assoluto (magnitudo)
    % 9. Fusione dei Detail Layers (Choose-Max basato su valore assoluto, Eq. 10)
    detailsA = coeffsA(:, :, 1:idx_base-1);
    detailsB = coeffsB(:, :, 1:idx_base-1);

    mask = abs(detailsA) >= abs(detailsB);
    coeffsF(:, :, 1:idx_base-1) = detailsA .* mask + detailsB .* (~mask);

    % (:, :, 1:idx_base-1) estrae in un colpo solo tutte le fette di
    % dettaglio (tutte tranne l'ultima, che è il base layer): con
    % nScales=4 e shearLevels=[1,1,2,2] sono le fette dalla 1 alla 48
    % (tutti i dettagli a tutte le 4 scale e direzioni). mask = ... confronta
    % il valore assoluto elemento per elemento. mask sarà una
    % matrice 3D booleana con valore 1 dove il dettaglio di A è più forte
    % di B e 0 altrove.

    % 10. Ricostruzione tramite Trasformata Shearlet Inversa (Sez. 3.4)
    imgF = SLshearrec2D(coeffsF, shearletSystem); % ricombina la base fusa e le 48 sotto-bande

    % Troncamento matematico nel range valido [0, 1]
    imgF = min(max(imgF, 0), 1);

    % 11. Salvataggio dell'immagine fusa su disco
    imwrite(imgF, fullfile(out_dir, ct_files(k).name));

    % Calcolo delle metriche di valutazione sulla coppia k-esima
    [Q, L, N, mad_val, sd_val] = compute_petrovic_metrics(imgA, imgB, imgF);
    metrics_results(k, :) = [Q, L, N, mad_val, sd_val];

end

fprintf('Elaborazione e fusione completata con successo per tutte le %d immagini!\n', num_images);


% =========================================================================
% STAMPA DELLA TABELLA DEI RISULTATI QUANTITATIVI
% =========================================================================
fprintf('\n=================================================================================\n');
fprintf('               TABELLA DEI RISULTATI QUANTITATIVI (Dataset CT-MRI)               \n');
fprintf('=================================================================================\n');
fprintf(' Coppia |   Q^{AB/F} (Info)  |   L^{AB/F} (Loss)  |   N^{AB/F} (Art)   |    MAD    |    SD     \n');
fprintf('---------------------------------------------------------------------------------\n');
for k = 1:num_images
    fprintf('  %02d    |      %0.4f       |      %0.4f       |      %0.4f       |  %0.4f  |  %0.4f   \n', ...
        k, metrics_results(k, 1), metrics_results(k, 2), metrics_results(k, 3), metrics_results(k, 4), metrics_results(k, 5));
end
fprintf('---------------------------------------------------------------------------------\n');
mean_metrics = mean(metrics_results, 1);
fprintf(' MEDIA  |      %0.4f       |      %0.4f       |      %0.4f       |  %0.4f  |  %0.4f   \n', ...
    mean_metrics(1), mean_metrics(2), mean_metrics(3), mean_metrics(4), mean_metrics(5));
fprintf('=================================================================================\n\n');