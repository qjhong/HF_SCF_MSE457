%normalization of basis set
function [coef_new]=normal(n_pg_basis,alpha_basis,coef_basis,type_basis)

    S=0;
    %Overlap integral
        for i=1:n_pg_basis
        for j=1:n_pg_basis
            S=S+overlap_int([0 0 0],alpha_basis(i),type_basis(:),[0 0 0],alpha_basis(j),type_basis(:))*coef_basis(i)*coef_basis(j);
        end
        end
        
        coef_new=coef_basis/sqrt(S);
end
    function [s]=overlap_int(center1,alpha1,type1,center2,alpha2,type2)
        P=(alpha1*center1+alpha2*center2)/(alpha1+alpha2);
        PA=P-center1;
        PB=P-center2;
        AB=center1-center2;
        AB2=AB*AB';
        gamma=alpha1+alpha2;
        l1=type1(1);
        l2=type2(1);
        i=floor((l1+l2)/2);
        fx=0;
        for ii=0:i
            fx=fx+f(2*ii,l1,l2,PA(1),PB(1))*double_factorial(2*ii-1)/(2*gamma)^ii;
        end
        m1=type1(2);
        m2=type2(2);
        j=floor((m1+m2)/2);
        fy=0;
        for jj=0:j
            fy=fy+f(2*jj,m1,m2,PA(2),PB(2))*double_factorial(2*jj-1)/(2*gamma)^jj;
        end
        n1=type1(3);
        n2=type2(3);
        k=floor((n1+n2)/2);
        fz=0;
        for kk=0:k
            fz=fz+f(2*kk,n1,n2,PA(3),PB(3))*double_factorial(2*kk-1)/(2*gamma)^kk;
        end
        s=(pi/gamma)^(3/2)*exp(-alpha1*alpha2*AB2/gamma)*fx*fy*fz;
    end
    function [df]=double_factorial(i)
        if (i<=1)
            df=1;
        else
            ii=i;
            df=1;
            while (ii>1)
                df=df*ii;
                ii=ii-2;
            end
        end
    end
    function [f]=f(j,l1,l2,a,b)
        f=0;
        j1=j-l2;
        if (j1<0)
            j1=0;
        end
        for i=j1:l1    %in the term x^j, j=i+(j-i), i is from the 1st, (j-i) is from the 2nd
            f=f+c(i,l1)*a^(l1-i)*c(j-i,l2)*b^(l2-j+i);
        end
    end
    function [c]=c(i,j)  %calculate C^i_j
        c=1;
        for l=1:i
            c=c*(j-l+1)/l;
        end            
    end
