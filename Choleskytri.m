
function x = Choleskytri(dp, ds, b)    
    N = length(dp);      
    L = zeros(N, 1);     % Diagonale de L
    Ls = zeros(N-1, 1);  % Sous-diagonale de L

    % Étape 1 : Factorisation de Cholesky
    L(1) = sqrt(dp(1));  % Premier élément de la diagonale de L
    for i = 1:N-1
        Ls(i) = ds(i) / L(i);          % Calcul de la sous-diagonale de L
        L(i+1) = sqrt(dp(i+1) - Ls(i)^2);  % Calcul de la diagonale de L
    end
    % disp('diagonale principale de L:')
    % disp(L)
    % disp('sous diagonale de L:')
    % disp(Ls)


    % Étape 2 : Résolution par substitution avant pour Ly = b
    y = zeros(N, 1);    % Vecteur y
    y(1) = b(1) / L(1);  % Premier élément de y
    for i = 2:N
        y(i) = (b(i) - Ls(i-1) * y(i-1)) / L(i);  % Résolution pour y
    end

    % Étape 3 : Résolution par substitution arrière pour L^T x = y
    x = zeros(N, 1);    % Vecteur x (solution du système)
    x(N) = y(N) / L(N); 
    for i = N-1:-1:1
        x(i) = (y(i) - Ls(i) * x(i+1)) / L(i);
    end
end
