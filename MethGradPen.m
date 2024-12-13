function [x, iter, err, p] = MethGradPen(A, b, x0, p, rho, tol, max_iter)
    % Algorithme du gradient avec pénalisation
    % 
    % Entrées:
    %   A: matrice du système
    %   b: second membre
    %   x0: solution initiale
    %   p: paramètre de pénalisation
    %   rho: paramètre de pénalisation
    %   tol: tolérance
    %   max_iter: nombre maximal d'itérations
    % Sorties:
    %   x: solution
    %   iter: nombre d'itérations
    %   tol: erreur relative
    %   p: pénalisation
    
    % Initialisation
    x = x0;  % Solution initiale
    iter = 0;  % Compteur d'itérations
    
    while iter < max_iter
        iter = iter + 1;
        
        % Calcul du gradient
        grad = A * x - b + 2 * p * norm(x) .* x;
        
        % Mise à jour de l'itération
        xt = x - rho * grad;
        
        % Calcul de l'erreur
        err = norm(xt - x);

        % Mise à jour de la solution
        x = xt;

        if err < tol
            break;
        end
    end
end