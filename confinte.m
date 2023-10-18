function [E_ci]=confinte(H,T,Va,S,V4,C,nzeta,fid)

    temp=size(C);
    n_basis=temp(1);
    %number ci matrix index
    nlist=1;
    list(:,nlist)=[1 1];    %[1 1] |C1aC1b>
    for i=2:nzeta           %S start
        nlist=nlist+1;
        list(:,nlist)=[1 i];
    end
    for i=2:nzeta
        nlist=nlist+1;
        list(:,nlist)=[i,1];
    end                     %S end
    for i=2:nzeta           %D start
    for j=2:nzeta
        nlist=nlist+1;
        list(:,nlist)=[i,j];
    end
    end                     %D end

    CIMATRIX=zeros(nlist,nlist);
    %evaluate matrix element of ci matrix
    H_CI=zeros(nzeta,nzeta);
    V4_CI=zeros(nzeta,nzeta,nzeta,nzeta);
    for i=1:nzeta
    for j=i:nzeta
        H_CI(i,j)=C(:,i)'*H*C(:,j);
        H_CI(j,i)=H_CI(i,j);
        T_CI(i,j)=C(:,i)'*T*C(:,j);
        T_CI(j,i)=T_CI(i,j);
        Va_CI(i,j)=C(:,i)'*Va*C(:,j);
        Va_CI(j,i)=Va_CI(i,j);
    end
    end
    %Psi1=|i1j1> Psi2=|i2j2>
    %<Psi1||Psi2>=(i1i2||j1j2)
    for i1=1:nzeta
    for j1=1:nzeta
    for i2=i1:nzeta
    for j2=j1:nzeta
        %V4(i,j,k,l): (ij||kl)
        %V4_CI(i1,j1,i2,j2): <i1j1||i2j2>=(i1i2||j1j2)
        if (abs(V4_CI(i1,j1,i2,j2))<1e-10)
            for m1=1:n_basis
            for m2=1:n_basis
            for n1=1:n_basis
            for n2=1:n_basis
                V4_CI(i1,j1,i2,j2)=V4_CI(i1,j1,i2,j2)+V4(m1,n1,m2,n2)*C(m1,i1)*C(n1,i2)*C(m2,j1)*C(n2,j2);
            end
            end
            end
            end
            V4_CI(i1,j2,i2,j1)=V4_CI(i1,j1,i2,j2);
            V4_CI(i2,j1,i1,j2)=V4_CI(i1,j1,i2,j2);
            V4_CI(i2,j2,i1,j1)=V4_CI(i1,j1,i2,j2);
            V4_CI(j1,i1,j2,i2)=V4_CI(i1,j1,i2,j2);
            V4_CI(j1,i2,j2,i1)=V4_CI(i1,j1,i2,j2);
            V4_CI(j2,i1,j1,i2)=V4_CI(i1,j1,i2,j2);
            V4_CI(j2,i2,j1,i1)=V4_CI(i1,j1,i2,j2);
        end
    end
    end
    end
    end
    
    for i=1:nlist
    for j=i:nlist
        i1=list(1,i);
        i2=list(2,i);
        i3=list(1,j);
        i4=list(2,j);
        %h1
        h1=0;
        if (i2==i4)
            h1=h1+H_CI(i1,i3);
        end
        if (i1==i3)
            h1=h1+H_CI(i2,i4);
        end
        h1=h1/2;
        %h2
        h2=0;
        if (i1==i3)
            h2=h2+H_CI(i2,i4);
        end
        if (i2==i4)
            h2=h2+H_CI(i1,i3);
        end
        h2=h2/2;
        V2=V4_CI(i1,i2,i3,i4);
        CIMATRIX(i,j)=h1+h2+V2;
        CIMATRIX(j,i)=CIMATRIX(i,j);
        CIMATRIX_Vee(i,j)=V2;
        CIMATRIX_Vee(j,i)=CIMATRIX_Vee(i,j);
        %h1
        h1=0;
        if (i2==i4)
            h1=h1+T_CI(i1,i3);
        end
        if (i1==i3)
            h1=h1+T_CI(i2,i4);
        end
        h1=h1/2;
        %h2
        h2=0;
        if (i1==i3)
            h2=h2+T_CI(i2,i4);
        end
        if (i2==i4)
            h2=h2+T_CI(i1,i3);
        end
        h2=h2/2;
        CIMATRIX_T(i,j)=h1+h2;
        CIMATRIX_T(j,i)=CIMATRIX_T(i,j);
        %h1
        h1=0;
        if (i2==i4)
            h1=h1+Va_CI(i1,i3);
        end
        if (i1==i3)
            h1=h1+Va_CI(i2,i4);
        end
        h1=h1/2;
        %h2
        h2=0;
        if (i1==i3)
            h2=h2+Va_CI(i2,i4);
        end
        if (i2==i4)
            h2=h2+Va_CI(i1,i3);
        end
        h2=h2/2;
        CIMATRIX_Va(i,j)=h1+h2;
        CIMATRIX_Va(j,i)=CIMATRIX_Va(i,j);
        %V2
%         V2=V4_CI(i1,i2,i3,i4);
%         for m1=1:nzeta
%         for m2=1:nzeta
%         for n1=1:nzeta
%         for n2=1:nzeta
%             V2=V2+V4(m1,n1,m2,n2)*C(m1,i1)*C(n1,i3)*C(m2,i2)*C(n2,i4);
%         end
%         end
%         end
%         end
%         if (abs(V2-V4_CI(i1,i2,i3,i4))>1e-5)
%             keyboard
%         end
    end
    end

    [C_CI,E_CI]=eig(CIMATRIX);
    for i=1:nlist
        for j=1:nlist-1
            if (E_CI(j+1,j+1)<E_CI(j,j))
                temp=E_CI(j,j);
                E_CI(j,j)=E_CI(j+1,j+1);
                E_CI(j+1,j+1)=temp;
                vtemp=C_CI(:,j);
                C_CI(:,j)=C_CI(:,j+1);
                C_CI(:,j+1)=vtemp;
            end
        end
    end
    
    CIMATRIX(1,1);
    E_CI(1:2,1:2);
    C_CI(:,1)';
    E_ci=C_CI(:,1)'*CIMATRIX(:,:)*C_CI(:,1);
    E_ci
    fprintf(fid,'%10.6f\n',E_CI(1,1));
    E_t=C_CI(:,1)'*CIMATRIX_T*C_CI(:,1);
    E_Va=C_CI(:,1)'*CIMATRIX_Va*C_CI(:,1);
    E_Vee=C_CI(:,1)'*CIMATRIX_Vee*C_CI(:,1);
    fprintf(fid,'CI kinetic energy: %10.6f\n',E_t);
    fprintf(fid,'CI external energy: %10.6f\n',E_Va);
    fprintf(fid,'CI internal energy: %10.6f\n',E_Vee);
    fprintf(fid,'CI total energy: %10.6f\n',E_ci);    
    