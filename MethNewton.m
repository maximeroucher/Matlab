function [x, iter] = MethNewton(A, b, x0, tol, max_iter)
    % Algorithme de Newton
    %
    % Entrées:
    %   A: matrice du système
    %   b: second membre
    %   x0: solution initiale
    %   tol: tolérance
    %   max_iter: nombre maximal d'itérations
    % Sorties:
    %   x: solution
    %   iter: nombre d'itérations

    % Initialisation
    x = x0;  % Solution initiale
    iter = 0;  % Compteur d'itérations
    d = length(b);  % Dimension du problème
    lambda = 0; % Initialisation de lambda
    u = [x; lambda]; % Initialisation de u

    
    while iter < max_iter
        jac = [A + 2 * lambda * eye(d), 2 * x; -2 * x', 0]; % Calcul de la jacobienne
        scmb = [A * x - b + 2 * lambda * x; -x' * x - 1]; % Calcul du second membre
        ut = u - jac\scmb; % iteration
        err = norm(u-ut); % calcul de l'erreur
        u = ut; % mise a jour
        iter = iter+1; % increment
        x=u(1:d,1); % extraction x
        lambda=u(d+1,1); % extraction lambda
        if err < tol % condition d'arret
            break;
        end
    end
end