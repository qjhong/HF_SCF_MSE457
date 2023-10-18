%6-311**G calculation on atom
%currently only for He, 31**G
%IOP=0  No printing whatsoever
%IOP=1  Print only converged results
%IOP=2  Print every iteration
function main

    %format of output
    IOP=2;
    
    %%%%%%%%% basic information here %%%%%%%%%%%%
    %atoms
    za=2;
    %basis set
    fidbasis=fopen('basis.txt','r');
    n_basis=fscanf(fidbasis,'%d',[1 1]);%total number of basis
    n_pg_basis=zeros(n_basis,1);    %number of premitive Gaussians in each basis
    alpha_basis=zeros(n_basis,6);   %exponent in Gaussian
    coef_basis=zeros(n_basis,6);    %coefficient in Gaussian
    center_basis=zeros(n_basis,3);       %center of each basis
    type_basis=zeros(n_basis,3);    %s,p,d,f
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FORMAT of basis.txt
    %1st: total number of basis
    %--------
    %1st: number of premitive gaussian in the basis
    %2nd: alpha of each gaussian
    %3rd: coef of each gaussian
    %4th: center of the basis
    %5th: type of the basis
    %--------
    for i=1:n_basis
        n_pg_basis(i)=fscanf(fidbasis,'%d',[1 1]);
        for j=1:n_pg_basis(i)
            alpha_basis(i,j)=fscanf(fidbasis,'%f',[1 1]);
        end
        for j=1:n_pg_basis(i)
            coef_basis(i,j)=fscanf(fidbasis,'%f',[1 1]);
        end
        for j=1:3
            center_basis(i,j)=fscanf(fidbasis,'%f',[1,1]);
        end
        for j=1:3
            type_basis(i,j)=fscanf(fidbasis,'%f',[1,1])';
        end
    end
    fclose(fidbasis);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Hartree-Fock calculation
    fid=fopen('result.txt','w');

    if (IOP~=0)
            fprintf(fid,'6-31**G for He\n');
    end
    %normalization of the basis set
    for i=1:n_basis
        coef_basis(i,:)=normal(n_pg_basis(i),alpha_basis(i,:),coef_basis(i,:),type_basis(i,:));
    end
    
    %calculate all the one and two-election integrals
    [S T Va V4]=intgrl(IOP,n_basis,n_pg_basis,alpha_basis,coef_basis,center_basis,type_basis,2,fid);
    %be inefficient and put all integrals in a pretty arrays
    [H,X]=colect(IOP,S,T,Va,fid);
    %perform the SCF calculation
    [C,ent]=SCF(IOP,H,T,Va,V4,X,n_basis,fid);
    %perform CI calculation
    confinte(H,T,Va,S,V4,C,n_basis,fid)

end