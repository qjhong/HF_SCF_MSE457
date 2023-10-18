function [C,ent]=SCF(IOP,H,T,Va,V4,X,n_basis,fid)
%perform the scf iterations

    %convergence criterion for density matrix
    crit=1.0e-4;
    %maximum number of iterations
    maxit=1000;
    %iteration number
    iter=0;
    stop=0;
    nsize=n_basis;
    
    %use core-hamiltonian for initial guess of F. I.E. (P=0)
    P=zeros(nsize,nsize);
%     P(3:4,3:4)=[  1.286143  0.540172;0.540172  0.226869];
    if (IOP==2)
        fprintf(fid,'P\n');
        for i=1:nsize
        for j=1:nsize
            fprintf(fid,'%10.6f',P(i,j));
        end
        fprintf(fid,'\n');
        end
    end
    %start of iteration loop
    while (stop==0)
        iter=iter+1;
        if (IOP==2)
            fprintf(fid,'start of iteration number %i\n',iter);
        end
        %form two-electron part of Fock matrix from P
        G=FORMG(P,V4);
        if (IOP==2)
            fprintf(fid,'G\n');
            for i=1:nsize
            for j=1:nsize
                fprintf(fid,'%10.6f',G(i,j));
            end
            fprintf(fid,'\n');
            end
        end
        %add core hamiltonian to get Fock matrix
        F=H+G;
        %calculate electronic energy
        en=0.0;
        en=en+0.5*sum(sum(P.*(H+F)));
        if (IOP==2)
            fprintf(fid,'F\n');
            for i=1:nsize
            for j=1:nsize
                fprintf(fid,'%10.6f',F(i,j));
            end
            fprintf(fid,'\n');
            end
            fprintf(fid,'electronic energy = %10.6f\n',en);
        end
        F1=X'*F*X;
        %diagonalize transformed Fock matrix
        [V D]=eig(F1);
        for i=1:nsize-1
        for j=1:nsize-i
            if (D(j,j)>D(j+1,j+1))
                temp=D(j,j);
                D(j,j)=D(j+1,j+1);
                D(j+1,j+1)=temp;
                tempv=V(:,j);
                V(:,j)=V(:,j+1);
                V(:,j+1)=tempv;
            end
        end
        end
        C=X*V;
        
        %form new density matrix
        %save present density matrix before creating a new one
        OLDP=P;
        for i=1:nsize
        for j=1:nsize
            P(i,j)=0.0;
            for k=1:1
                P(i,j)=P(i,j)+2.0*C(i,k)*C(j,k);
            end
        end
        end

        if (IOP==2)
            fprintf(fid,'EWs\n');
            fprintf(fid,'%10.6f%10.6f\n',D(1,1),D(2,2));
            fprintf(fid,'C\n');
            for i=1:nsize
            for j=1:nsize
                fprintf(fid,'%10.6f',C(i,j));
            end
            fprintf(fid,'\n');
            end
            fprintf(fid,'P\n');
            for i=1:nsize
            for j=1:nsize
                fprintf(fid,'%10.6f',P(i,j));
            end
            fprintf(fid,'\n');
            end
            Ek=C(:,1)'*T*C(:,1)*2;
            Vext=C(:,1)'*Va*C(:,1)*2;
            Vhartree=0.0;
            for i=1:nsize
            for j=1:nsize
            for k=1:nsize
            for l=1:nsize
                Vhartree=Vhartree+C(i,1)*C(j,1)*C(k,1)*C(l,1)*V4(i,j,k,l);
            end
            end
            end
            end
            Vhartree=Vhartree/2*4;
            fprintf(fid,'\nkinetic energy is %10.6f',Ek);
            fprintf(fid,'\nexternal energy is %10.6f',Vext);
            fprintf(fid,'\nhartree energy is %10.6f',Vhartree);
            fprintf(fid,'\nexchange energy is %10.6f',en-Ek-Vext-Vhartree);
            fprintf(fid,'\ntotal energy is %10.6f\n',en);
        end        
        P=0.9*OLDP+0.1*P;
        if (IOP==2)
            fprintf(fid,'EWs\n');
            fprintf(fid,'%10.6f%10.6f\n',D(1,1),D(2,2));
            fprintf(fid,'C\n');
            for i=1:nsize
            for j=1:nsize
                fprintf(fid,'%10.6f',C(i,j));
            end
            fprintf(fid,'\n');
            end
            fprintf(fid,'P\n');
            for i=1:nsize
            for j=1:nsize
                fprintf(fid,'%10.6f',P(i,j));
            end
            fprintf(fid,'\n');
            end
        end
        %calculate delta
        delta=sum(sum((P-OLDP).^2));
        delta=sqrt(delta/4);
        if (IOP>0)
            fprintf(fid,'convergence of density matirx = %10.6f\n',delta);
        end
        %check for convergence
        if (delta<crit)
            %calculation converged if it got here
            %add nuclear repulsion to get total energy
%             ent=en+za*zb/r;
            ent=en;
            if (IOP>0)
                fprintf(fid,'calculation converged. electronic energy = %10.6f, total energy = %10.6f\n', en, ent);
            end
            %print out the final results
            if (IOP==2)
                fprintf(fid,'2%10.6f',G);
                fprintf(fid,'2%10.6f',F);
%                 fprintf(fid,'2%10.6f',E)
                fprintf(fid,'2%10.6f',C);
                fprintf(fid,'2%10.6f',P);
            end

            stop=1;
        end
        %not yet converged
        %test for maximum number of iterations
        %if maximum number not yet reached then go back for another iteration
        if (iter>maxit)
            fprintf(fid,'something wrong here');
            stop=1;
        end
    end
        
end
function [G]=FORMG(P,V4)
%calculates the G matrix from the density matrix and two-electron integrals
    temp=size(P);
    nsize=temp(1);
    G=zeros(nsize,nsize);
    for i=1:nsize
    for j=1:nsize
    for k=1:nsize
    for l=1:nsize
        G(i,j)=G(i,j)+P(k,l)*(V4(i,j,k,l)-0.5D0*V4(i,l,k,j));
    end
    end
    end
    end
end