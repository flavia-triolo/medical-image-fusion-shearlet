function q = guided_filter(I, p, r, eps_val)
% 1. Definizione della finestra locale quadrata (2r+1) x (2r+1)
w_size = 2 * r + 1;
h = ones(w_size, w_size) / (w_size^2);
% se il raggio è r=30, la larghezza
% della finestra locale è w_size = 2*r+1=61 pixel. Il comando ones(61, 61)
% crea una matrice 61x61 in cui tutti gli elementi valgono 1. dividendo per
% l'area totale della finestra (w_size^2=61^2=3721), si ottiene la matrice
% h: maschera di media
% locale. Poichè la somma di tutti i suoi coefficienti fa esattamente 1,
% fare una convoluzione con h equivale a fare la media aritmetica classica
% di tutti i pixel dentro quel quadratino 61x61
% 2. Calcolo delle medie locali E[I], E[p], E[I*p]
mean_I  = imfilter(I, h, 'replicate');      %è la media locale del''immagine guida che coincide con il valore atteso E(X) rispetto alla misura di probabilitò uniforme definita sulla finestra w
mean_p  = imfilter(p, h, 'replicate');      % E[p]
mean_Ip = imfilter(I .* p, h, 'replicate'); %E[I*p]

%imfilter è la funzione MATLAB che calcola la convoluzione spaziale 2D tra
%una matrice X e il filtro h. Prende la finestra h di 61x61pixel e la
%centra sul pixel (x,y) dell'immagine X. moltiplica ogni pixel del vicinato
%e somma tutto. L'opzione 'replicate' serve ai bordi esterni dell'immagine:
%quando la finestra 61x61 "sborda" fuori dall'immagine, MATLAB replica i
%valori dei bordi anzichè mettere zeri neri artificiali.

% 3. Covarianza locale Cov(I, p) = E[I*p] - E[I]*E[p]
cov_Ip = mean_Ip - mean_I .* mean_p;  %Cov(X,Y)=E[XY]-E[X]E[Y]. misura quanto la mappa binaria varia concordamente con la guida

% 4. Varianza locale Var(I) = E[I^2] - (E[I])^2
mean_II = imfilter(I .* I, h, 'replicate');   %E[I^2]
var_I   = mean_II - mean_I .* mean_I;   %calcola la varianza locale

% 5. Coefficienti di regressione locale a e b (soluzione analitica del minimo)
a = cov_Ip ./ (var_I + eps_val);
b = mean_p - a .* mean_I;

% 6. Media dei coefficienti a e b sulle finestre sovrapposte
mean_a = imfilter(a, h, 'replicate');
mean_b = imfilter(b, h, 'replicate');

% 7. Output finale filtrato: q = mean_a * I + mean_b
q = mean_a .* I + mean_b;
end
