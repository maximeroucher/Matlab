function [x, y, J_values, norm_values, iter] = MethGradPO(max_iter)
    % ======================================================================
    % Cette fonction implémente la méthode du gradient à pas optimal pour
    % minimiser la fonction J.
    %
    % Paramètres :
    %   max_iter : Nombre maximal d'itérations
    %
    % Sortie :
    %   x : Coordonnée x du minimum de la fonction J
    %   y : Coordonnée y du minimum de la fonction J
    %   J_values : Valeurs de la fonction J à chaque itération
    %   norm_values : Normes de la différence Xk - X* à chaque itération
    %   iter : Nombre d'itérations effectuées
    % ======================================================================
    % Initialisation des variables
    x = 0; % Point initial en (x0, y0) = (0, 0)
    y = 0; 
    J_values = []; % Stocker les valeurs de la fonction J
    norm_values = []; % Stocker les normes de la différence Xk - X*
    iter = 0; % Initialisation du compteur d'itérations
    X_star = [1, 1]; % Point de minimum de la fonction J
    
    % Boucle d'optimisation (itération du gradient)
    while iter < max_iter
         iter = iter + 1;
        
        % Calcul du gradient à (x, y)
        [grad_x, grad_y] = GradJ([x, y]);
        
        % Calcul du pas optimal avec la fonction pasOptimal
        r0 = 0.1; % Point initial pour la recherche du pas optimal
        X = [x, y]; % Point actuel (x, y)
        
        % Appel à la fonction pasOptimal pour obtenir alpha optimal et fval
        [alpha_opt, fval] = pasOptimal(r0, X);
        
        % Mise à jour des variables x et y en fonction du pas optimal
        x = x - alpha_opt * grad_x;
        y = y - alpha_opt * grad_y;
        
        % Calcul de la valeur de la fonction J à ce point
        J_val = (x - 1)^2 + 10 * (x^2 - y)^2;
        J_values = [J_values, J_val]; % Stockage de la valeur de J
        
        % Calcul de la norme de Xk - X*
        norm_value = norm([x - X_star(1), y - X_star(2)]); % Calcul de la norme
        norm_values = [norm_values, log(norm_value)]; % Stockage du log de la norme
    end
end