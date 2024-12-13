function [x, iter] = MethGradProj(A, b, x0, rho, tol, max_iter)
    % Algorithme du gradient avec projection sur la sphère unité
    %
    % Entrées:
    %   A: matrice du système
    %   b: second membre
    %   x0: solution initiale
    %   rho: pas de gradient
    %   tol: tolérance
    %   max_iter: nombre maximal d'itérations
    % Sorties:
    %   x: solution
    %   iter: nombre d'itérations

    
    % Initialisation
    x = x0;  % Solution initiale
    iter = 0;  % Compteur d'itérations
    d = length(b);  % Dimension du problème
    
    while iter < max_iter
        iter = iter + 1;
        
        % Calcul du gradient
        grad = A * x - b;
        
        % Mise à jour de l'itération
        x_new = x - rho * grad;
        
        % Projection sur la sphère unité
        if norm(x_new) > 1
            x = x_new / norm(x_new);  % Normalisation si nécessaire
        else
            x = x_new;  % Sinon, on conserve la solution
        end

        % Vérification de la condition d'arrêt
        if norm(x - x_new) < tol
            break;
        end
    end
end