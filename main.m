%%  A
% Appel de la fonction pour afficher le graphe sur le domaine [-10, 10]
a = 10;    % Limites de la zone [-10, 10]
h = 0.1;   % Pas de discrétisation
% plotGraphe(a, h);

% Appel de la fonction pour afficher les lignes de niveau sur le domaine [-10, 10]
a = 10;    % Limites de la zone [-10, 10]
d = 0.1;   % Pas de discrétisation
N = 10;    % Nombre de lignes de niveau
% LignesNiveau(a, d, N);

%% B - I

% Tracés
max_iter = 200; % Nombre maximal d'itérations
alpha = 0.01; % Pas constant pour la méthode du gradient

% Méthode du Gradient à Pas Constant
[x1, y1, J_values1, norm_values1, iter1] = MethGrad(alpha, max_iter);

% Méthode du Gradient à Pas Optimal
[x2, y2, J_values2, norm_values2, iter2] = MethGradPO(max_iter);

% Tracé des courbes de convergence
figure;
plot(1:iter1, norm_values1, 'r', 'LineWidth', 2); % Méthode du gradient à pas constant
hold on;
plot(1:iter2, norm_values2, 'b', 'LineWidth', 2); % Méthode du gradient à pas optimal
xlabel('Itérations k');
ylabel('log(\|X_k - X^*\|)');
legend('Gradient à Pas Constant', 'Gradient à Pas Optimal');
title('Convergence des Méthodes du Gradient');
grid on;

%% B - II


function [x, iter, errors] = GradOpt(A, b, x0, tol, max_iter)
    % Algorithme du gradient à pas optimal avec suivi de l'erreur relative
    
    % Initialisation
    x = x0;  
    r = b - A * x;  
    iter = 0;  
    r0_norm = norm(r, inf);  % Calcul de la norme infinie du résidu initial
    errors = [];  % Liste pour stocker les erreurs relatives à chaque itération
    
    while iter < max_iter
        iter = iter + 1;
        
        % Calcul du gradient
        grad = -r;  % Le gradient de J est simplement le résidu
        
        % Calcul du pas optimal (alpha optimal)
        phi = @(alpha) norm(b - A * (x + alpha * grad), 2)^2;  % Fonction à minimiser
        [alpha_optimal, fval] = fminsearch(phi, 0);  % Recherche du minimum de la fonction
        
        % Mise à jour de la solution
        x = x + alpha_optimal * grad;
        
        % Mise à jour du résidu
        r_new = b - A * x;
        
        % Calcul de l'erreur relative
        error = norm(r_new, inf) / r0_norm;
        errors = [errors, error];  % Stocker l'erreur relative
        
        % Vérification de la condition d'arrêt
        if norm(r_new, inf) < tol
            break;
        end
        
        % Mise à jour du résidu
        r = r_new;
    end
end

% Paramètres
n = 20;
A = matdiff(n);
x0 = zeros(n, 1);  % Solution initiale
tol_conj = 1e-6;  % Tolérance pour le critère d'arrêt du gradient conjugué
tol_opt = 1e-3;  % Tolérance pour le critère d'arrêt du gradient à pas optimal
max_iter = 200;  % Nombre maximal d'itérations

% Cas 1: Second membre b1
b1 = [0; ones(n-2, 1); 0];

% Gradient Conjugué pour b1
[x1_gc, iter1_gc, errors1_gc] = GradConj(A, b1, x0, tol_conj, max_iter);

% Gradient à pas optimal pour b1
[x1_opt, iter1_opt, errors1_opt] = GradOpt(A, b1, x0, tol_opt, max_iter);

% Cas 2: Second membre b2
b2 = rand(n, 1);

% Gradient Conjugué pour b2
[x2_gc, iter2_gc, errors2_gc] = GradConj(A, b2, x0, tol_conj, max_iter);

% Gradient à pas optimal pour b2
[x2_opt, iter2_opt, errors2_opt] = GradOpt(A, b2, x0, tol_opt, max_iter);

% Comparaison des erreurs pour le cas b1 et b2
figure;
semilogy(1:length(errors1_gc), errors1_gc, 'r-', 'LineWidth', 2); hold on;
semilogy(1:length(errors1_opt), errors1_opt, 'r--', 'LineWidth', 2);
semilogy(1:length(errors2_gc), errors2_gc, 'b-', 'LineWidth', 2);
semilogy(1:length(errors2_opt), errors2_opt, 'b--', 'LineWidth', 2);

xlabel('Itérations');
ylabel('Erreur relative (log)');
legend('Gradient conjugué b1', 'Gradient à pas optimal b1', 'Gradient conjugué b2', 'Gradient à pas optimal b2');
title('Comparaison des erreurs relatives pour le gradient conjugué et le gradient à pas optimal');
grid on;

%% B - III

% Paramètres
n = 20;  % Taille de la matrice
A = hilb(n);  % Matrice de Hilbert
x0 = zeros(n, 1);  % Solution initiale
b = ones(n, 1);  % Second membre
tol = 1e-12;  % Tolérance
max_iter = 200;  % Nombre maximal d'itérations

% Gradient Conjugué pour Hilbert
[x_gc, iter_gc, errors_gc] = GradConj(A, b, x0, tol, max_iter);

% Gradient à pas optimal pour Hilbert
[x_opt, iter_opt, errors_opt] = GradOpt(A, b, x0, tol, max_iter);

% Tracer des courbes semi-logarithmiques de l'erreur relative
figure;
semilogy(1:iter_gc, errors_gc, 'r-', 'LineWidth', 2); hold on;
semilogy(1:iter_opt, errors_opt, 'b-', 'LineWidth', 2);

xlabel('Itérations');
ylabel('Erreur relative (log)');
legend('Gradient conjugué', 'Gradient à pas optimal');
title('Comparaison des erreurs relatives pour le gradient conjugué et le gradient à pas optimal avec la matrice de Hilbert');
grid on;

%% C - I

% Paramètres
d_values = [10, 20, 50];  % Différentes tailles de la matrice A
max_iter = 200;  % Nombre maximal d'itérations
tol = 1e-6;  % Tolérance pour le critère d'arrêt

% Test de la méthode pour différentes tailles de la matrice A
for d = d_values
    % Création de la matrice de diffusion A (tridiagonale)
    A = matdiff(d);
    b = ones(d, 1);  % Second membre b = ones(d, 1)
    x0 = zeros(d, 1);  % Point initial (solution initiale)
    % Test pour différentes valeurs de rho
    rho_values = [0.1, 0.05, 0.01, 0.005];

    for rho = rho_values
        [x, iter] = MethGradProj(A, b, x0, rho, tol, max_iter);
        fprintf('Pour rho = %.4f, le nombre d''itérations est : %d\n', rho, iter);
        % Affichage du résultat
        fprintf('Pour d = %d, le nombre d''itérations est : %d\n', d, iter);
        
        % Calcul de l'erreur finale (norme de la différence avec la solution exacte)
        x_exact = A \ b;  % Solution exacte du problème linéaire
        error = norm(x - x_exact);
        fprintf('Erreur pour d = %d : %.4e\n', d, error);
    end
end

%% C - II

% Paramètres
d = 5; % Taille de la matrice A
A = matdiff(d);  % Matrice de diffusion
rho = [0.1, 0.05, 0.01, 0.005];  % Valeurs de rho
tol = 1e-6;  % Tolérance pour le critère d'arrêt
p = [1, 10, 100];  % Valeurs de p
b = ones(d, 1);  % Second membre

% Test de la méthode pour différentes valeurs de p et rho
for p_value = p
    for rho_value = rho
        % Initialisation de la solution x0
        x0 = zeros(d, 1);
        % Appel de la fonction MethGradPen
        [x, iter, err, p] = MethGradPen(A, b, x0, p, rho_value, tol, max_iter);
        % Affichage des résultats
        fprintf('Pour p = %d et rho = %.4f, le nombre d''itérations est : %d\n', p_value, rho_value, iter);
        % Calcul de l'erreur finale (norme de la différence avec la solution exacte)
        x_exact = A \ b;  % Solution exacte du problème linéaire
        error = norm(x - x_exact);
        fprintf('Erreur pour p = %d : %.4e\n', p_value, error);
    end
end

%% C - III

% Paramètres
d = 5; % Taille de la matrice A
A = matdiff(d);  % Matrice de diffusion
rho = [0.1, 0.05, 0.01, 0.005];  % Valeurs de rho
tol = 1e-6;  % Tolérance pour le critère d'arrêt

% Test de la méthode pour différentes valeurs de rho
for rho_value = rho
    % Initialisation de la solution x0
    x0 = zeros(d, 1);
    % Appel de la fonction MethNewton
    [x, iter] = MethNewton(A, b, x0, tol, max_iter);
    % Affichage des résultats
    fprintf('Pour rho = %.4f, le nombre d''itérations est : %d\n', rho_value, iter);
    % Calcul de l'erreur finale (norme de la différence avec la solution exacte)
    x_exact = A \ b;  % Solution exacte du problème linéaire
    error = norm(x - x_exact);
    fprintf('Erreur pour rho = %.4f : %.4e\n', rho_value, error);
end