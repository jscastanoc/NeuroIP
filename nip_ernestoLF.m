function [L, cortex] = nip_ernestoLF(lpath,head)

Lx = load(strcat(lpath,head,'_Pots_X.txt'));
Ly = load(strcat(lpath,head,'_Pots_Y.txt'));
Lz = load(strcat(lpath,head,'_Pots_Z.txt'));
L = cat(3,Lx',Ly',Lz');
L = nip_translf(L);

fname = strcat(lpath,head,'_vc.ply');
cortex.vc = importdata(fname);

fname = strcat(lpath,head,'_tri.ply');
tri = importdata(fname);
cortex.tri = tri(:,2:end)+1;


end