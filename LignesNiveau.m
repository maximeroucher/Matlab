function [] = LignesNiveau(a,d,N)
%______________________________________________________________
% a : borne du pave
% delta :pas de discretisation
% N : nombre de lignes de niveau
%______________________________________________________________
%
x=-a:d:a;
y=-a:d:a;
%
[X,Y]=meshgrid(x,y);
%
% Exple 1 : J(x,y) = 2xy - 4x + x^2 + 2y^2
% ----------------------------------------
%Z = -2*X.*Y - 4*X + X.^2 + 2*Y.^2;
%
% Exple 2 : J(x,y) = (x-1)^2 + 10*(x^2 -y)^2
% ------------------------------------------
Z = (X -1).^2 + 10*(X.^2 - Y).^2;
%
[C,h]=contour(X, Y, Z, N);
set(h,'ShowText','off','TextStep',get(h,'LevelStep')*2)
colormap cool