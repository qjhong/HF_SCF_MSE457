%STO-3G calculation on atom
%currently only for Be, STO-3G
%IOP=0  No printing whatsoever
%IOP=1  Print only converged results
%IOP=2  Print every iteration
function main2
%format of output
    IOP=2;
    
    %%%%%%%%% basic information here %%%%%%%%%%%%
    %atom
    Za=4;
    %basis set
    fidbasis=fopen('basis_Be.txt','r');
    n_basis=fscanf(fidbasis,'%d',[1 1]); %total number of basis
    n_pg_basis=zeros(n_basis,1);    %number of premitive Gaussians in each basis
    alpha_basis=zeros(n_basis,6);   %exponent in Gaussian
    coef_basis=zeros(n_basis,6);    %coefficient in Gaussian
    center_basis=zeros(n_basis,3);  %center of each basis
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
    n_basis1s = n_basis;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Hartree-Fock calculation
    fid=fopen('result.txt','w');

    if (IOP~=0)
            fprintf(fid,'STO-3G for Be\n');
    end
    %normalization of the basis set
    for i=1:n_basis
        coef_basis(i,:)=normal(n_pg_basis(i),alpha_basis(i,:),coef_basis(i,:),type_basis(i,:));
    end
    
    %calculate all the one and two-election integrals
    [S, T, Va, V4]=intgrl1s(IOP,n_basis,n_pg_basis,alpha_basis,coef_basis,center_basis,type_basis,Za,fid);
    %be inefficient and put all integrals in a pretty arrays
    [H,X]=colect(IOP,S,T,Va,fid);
    %perform the SCF calculation
    [C1s,G1s, F1s, P1s, ent1s]=SCF1s(IOP,H,T,Va,V4,X,n_basis,fid);
    %perform CI calculation
    E_ci1s = confinte1s(H,T,Va,S,V4,C1s,n_basis,fid);
    %% 
    Za=2; %The 2s Shell
    fidbasis=fopen('basis_6-31G.txt','r');
    n_basis=fscanf(fidbasis,'%d',[1 1]); %total number of basis
    n_pg_basis=zeros(n_basis,1);    %number of premitive Gaussians in each basis
    alpha_basis=zeros(n_basis,6);   %exponent in Gaussian
    coef_basis=zeros(n_basis,6);    %coefficient in Gaussian
    center_basis=zeros(n_basis,3);  %center of each basis
    type_basis=zeros(n_basis,3);    %s,p,d,f
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
    n_basis2s = n_basis;
    %Hartree-Fock calculation
%     fid=fopen('result_2s.txt','w');
% 
%     if (IOP~=0)
%             fprintf(fid,'STO-3G for Be\n');
%     end
    %normalization of the basis set
    for i=1:n_basis
        coef_basis(i,:)=normal(n_pg_basis(i),alpha_basis(i,:),coef_basis(i,:),type_basis(i,:));
    end
    
    %calculate all the one and two-election integrals
    [S, T, Va, V4]=intgrl2s(IOP,n_basis,n_pg_basis,alpha_basis,coef_basis,center_basis,type_basis,Za,fid);
    %be inefficient and put all integrals in a pretty arrays
    [H,X]=colect(IOP,S,T,Va,fid);
    %perform the SCF calculation
    [C2s, G2s, F2s, P2s, ent2s]=SCF2s(IOP,H,T,Va,V4,X,n_basis,fid);
    %perform CI calculation
    E_ci2s = confinte2s(H,T,Va,S,V4,C2s,n_basis,fid);
    HF_Be = E_ci1s+E_ci2s
    
    %Print final results
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'Final Results:\n');
    fprintf(fid,'G\n');
           for i=1:n_basis1s
           for j=1:n_basis1s
               if i>n_basis2s || j>n_basis2s % Note: This section only works for n_basis1s>n_basis2s
                   fprintf(fid,'%10.6f',G1s(i,j));
               else
                   fprintf(fid,'%10.6f',G1s(i,j)+G2s(i,j));
               end
           end
           fprintf(fid,'\n');
           end
    fprintf(fid,'F\n');
           for i=1:n_basis1s
           for j=1:n_basis1s
               if i>n_basis2s || j>n_basis2s
                   fprintf(fid,'%10.6f',F1s(i,j));
               else
                   fprintf(fid,'%10.6f',F1s(i,j)+F2s(i,j));
               end
           end
           fprintf(fid,'\n');
           end
    fprintf(fid,'C\n');
           for i=1:n_basis1s
           for j=1:n_basis1s
               if i>n_basis2s || j>n_basis2s
                   fprintf(fid,'%10.6f',C1s(i,j));
               else
                   fprintf(fid,'%10.6f',C1s(i,j)+C2s(i,j));
               end
           end
           fprintf(fid,'\n');
           end
    fprintf(fid,'P\n');
           for i=1:n_basis1s
           for j=1:n_basis1s
               if i>n_basis2s || j>n_basis2s
                   fprintf(fid,'%10.6f',P1s(i,j));
               else
                   fprintf(fid,'%10.6f',P1s(i,j)+P2s(i,j));
               end
           end
           fprintf(fid,'\n');
           end
    fprintf(fid,'total energy: %10.6f\n',ent1s + ent2s);
    fprintf(fid,'CI total energy: %10.6f\n',E_ci1s + E_ci2s);
end