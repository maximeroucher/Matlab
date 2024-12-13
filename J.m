function z = J(x, y)
    % Calcul de la fonction de Rosenbrock
    z = (x - 1).^2 + 10*(x.^2 - y).^2;
end