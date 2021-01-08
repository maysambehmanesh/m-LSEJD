function [V,L_sym,E] = laplacian_eigs_sym(W,eigsNum)

D = diag(sum(W)); 
L = D - W; 
normal=sparse(full(D)^(-1/2));
L_sym = normal * L * normal;

% make sure L_sym is symmetric ( not just numerically...)
L_sym = (L_sym+L_sym')/2;

% if you have probs with the shift value swap the two commands
[Vu,E] = eigs(L_sym,eigsNum,'SM');
%[Vu,E] = eigs(L_sym,eigsNum,1e-8);

% normalize V 

% NOTE: WHAT WE OBTAIN THIS WAY IS THE EIGENSPACE OF THE "RANDOM WALK" 
%       NORMALIZED LAPLACIAN -- i.e. (D^-1) * L
% V = normal * Vu;

% NOTE: WHAT WE OBTAIN THIS WAY IS THE (NORMALIZED) EIGENSPACE OF THE
%       "L_sym" SYMMETRIC LAPLACIAN
% V = bsxfun(@times, Vu, 1./(sum(Vu.^2, 2)));

% just return Vunnormalized
V = Vu;

[~,i] = sort(diag(E));
% return V ordered according to eigenvalues
V = V(:,i);
end