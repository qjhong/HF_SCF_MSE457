function [H,X]=colect(IOP,S,T,Va,fid)
%this takes the basic integrals from common and assembles the relevant
%matrices, that is S, H, X, XT, and two-electron integrals.

    %form core hamiltonian
    H=T+Va;
    %S is already formed.
    %form X: use canonical orthogonalization
    [V D]=eig(S);
    X=V*D^(-0.5)*V';
%     X(1,1)=1/sqrt(2*(1+S(1,2)));
%     X(2,1)=X(1,1);
%     X(1,2)=1/sqrt(2*(1-S(1,2)));
%     X(2,2)=-X(1,2);
    if (IOP~=0)
        temp=size(X);
        n_basis=temp(1);
        fprintf(fid,'H\n');
        for i=1:n_basis
        for j=1:n_basis
            fprintf(fid,'%10.6f',H(i,j));
        end
        fprintf(fid,'\n');
        end
        fprintf(fid,'X\n');
        for i=1:n_basis
        for j=1:n_basis
            fprintf(fid,'%10.6f',X(i,j));
        end
        fprintf(fid,'\n');
        end
    end

end