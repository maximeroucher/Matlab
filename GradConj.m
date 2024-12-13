function [x, iter, errors] = GradConj(A, b, x0, tol, max_iter)
    % ======================================================================
    % Cette fonction implémente la méthode du gradient conjugué pour résoudre
    % un système linéaire Ax = b.
    % 
    % Paramètres :
    %   A : Matrice du système linéaire
    %   b : Vecteur du second membre
    %   x0 : Point initial
    %   tol : Tolérance pour le critère d'arrêt
    %   max_iter : Nombre maximal d'itérations
    %
    % Sortie :
    %   x : Solution du système linéaire
    %   iter : Nombre d'itérations effectuées
    % ======================================================================
    
    % Initialisation
    x = x0;  % Point initial
    r = b - A * x;  % Résidu initial
    p = r;  % Direction initiale
    iter = 0;  % Compteur d'itérations
    r0_norm = norm(r, inf);  % Calcul de la norme infinie du résidu initial
    errors = [];  % Liste pour stocker les erreurs relatives à chaque itération
    
    while iter < max_iter
        iter = iter + 1;
        
        % Calcul du pas alpha
        alpha = (r' * r) / (p' * A * p);
        
        % Mise à jour de la solution
        x = x + alpha * p;
        
        % Mise à jour du résidu
        r_new = r - alpha * A * p;
        
        % Vérification de la condition d'arrêt
        if norm(r_new, inf) < tol
            break;
        end
        
        % Calcul du coefficient beta
        beta = (r_new' * r_new) / (r' * r);
        
        % Mise à jour de la direction
        p = r_new + beta * p;
        
        % Mise à jour du résidu
        r = r_new;

        % Calcul de l'erreur relative
        error = norm(r_new, inf) / r0_norm;
        errors = [errors, error];  % Stocker l'erreur relative
    end
end