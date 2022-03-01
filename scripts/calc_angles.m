function [alpha, beta, gamma] = calc_angles(coord_prox, coord_dist)

% elementi della matrice di rotazione (matrice (10) slide 10) utili per 
% calcolare gli angoli:
% - r31 = assez_prossimale * assex_distale % alpha -> rotazione caviglia-piede
% - r32 = assez_prossimale * assey_distale % beta
% - r12 = assex_prossimale * assey_distale
% elementi calcolati come prodotto scalare/norma

r31 = dot(coord_prox.z, coord_dist.x,2)./(vecnorm(coord_prox.z,2,2).*vecnorm(coord_dist.x,2,2));
r32 = dot(coord_prox.z, coord_dist.y,2)./(vecnorm(coord_prox.z,2,2).*vecnorm(coord_dist.y,2,2));
r12 = dot(coord_prox.x, coord_dist.y,2)./(vecnorm(coord_prox.x,2,2).*vecnorm(coord_dist.y,2,2));

% dopo aver calcolato gli angoli, trovare alpha, beta e gamma come segue:
% - alpha = arcsin(r_32)
% - beta = arcsin(-r_31/cos(alpha))
% - gamma = arcsin(-r_12/cos(alpha))
% alpha, beta e gamma sono rispettivamente l'angolo di:
% - adduzione-abduzione
% - flessione ed estensione
% - intra- ed extra-rotazione
% dovranno poi essere plottati singolarmente rispetto al loro piano di
% riferimento.

alpha = asin(r32); 
beta = asin(-r31./cos(alpha));
gamma = asin(-r12./cos(alpha));

end