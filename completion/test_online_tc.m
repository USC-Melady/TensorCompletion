sz = [5,8,9];
X_mat = rand(sz);
X = tensor(X_mat);


% completable subtensos
num_idx = 2;
omega_1 = [1, 4];
omega_2 = [2, 6];

omegas = cell(num_idx, 1);
omegas{1} = omega_1;
omegas{2} = omega_2;

% solve nuclear optimization
