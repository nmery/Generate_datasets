function [x0s,somme,poids,idx] = invd_new_3d(x,x0,a,dmax,nmax)
%
% Effectue l'estimation par inverse de la distance
%
% Synthaxe:	 [x0s,poids] = invd(x,x0,a,dmax)
%
%	x: matrice [x,y,z=valeur]contenant les données servant à l'estimation
%	x0: matrice contenant les coordonnées des point à estimer (n x 2)
%	a: exposant dans l'inverse de la distance
%   dmax : distance maximale
%   nmax : nombre de points maximal
%   dmax: rayon de recherche autour du point à estimer
%       x0s est x0 + valeur estimée en chaque point  
%       poids: la matrice des poids
%       idx: indice des points voisins
% 
nx=size(x,1);
nx0=size(x0,1);
nmax=min(nmax,nx);

[idx,dist]=knnsearch(x(:,1:3),x0,'K',nmax);

id=dist>dmax;
poids=zeros(nx0,nmax);

dist=dist.^a;
dist=1./dist;
dist(id)=0;

somme=sum(dist,2);
poids=dist./repmat(somme,1,nmax);
%clear somme dist 
est=zeros(nx0,1);
y=x(:,4);
for i=1:nx0;
   est(i)=poids(i,:)*y(idx(i,:));
end
x0s=[x0,est];

