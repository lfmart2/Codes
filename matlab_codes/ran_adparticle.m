clear variables
tic
N = 6000; v = 10; w = 10;
e0 = 0.0; T1 = 0.1; dlt = 2*0.05;
NN = 50; xi = 7;
nmol = linspace(1,100,100);
psol = 10;
el_oc = zeros([length(nmol) NN]);
for j = 1:length(nmol)
    for ii = 1:NN
        [H, mol] = solSSH_ran(N, v, w, e0, dlt, T1, xi, psol, nmol);
        [V, D] = eig(H);
        n00 = 0;
        for n = 1:length(mol)
            % Replace the inner loop with a matrix operation
            n00 = n00 + sum(V(N+n, 1:N/2).^2, 'all');
        end
        el_oc(j,ii) = n00;
    end
    j
end
noc = el_oc ./ (nmol');
save ran_part_psol10.mat noc w v N e0 T1 dlt NN xi nmol psol
toc