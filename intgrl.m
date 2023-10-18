%subroutine intgrl
%calculates all the basic integrals needed for SCF calculation
%parameters:
%n_basis(1): total number of basis
%n_pg_basis(n_basis,1): for each basis, the number of premitive Gaussians
%alpha_basis(n_basis,6): the exponential parameter in Gaussian
%coef_basis(n_basis,6): in each basis, coefficents of the premitive Gaussians
%center_basis(n_basis,3): the center of each basis
%type_basis(n_basis,3): x^i * y^j * z^k, i,j,k
function [S T Va V4]=intgrl(IOP,n_basis,n_pg_basis,alpha_basis,coef_basis,center_basis,type_basis,Za,fid)

    %read in the Ryspolynomial table
    fid2=fopen('Ryspoly.txt','r');
    Ryspolytable=zeros(1001,5,10);
    for i=1:5
    for j=1:1001
        x=0.1*(j-1);
        n=i;
        xx=fscanf(fid2,'%f %f',[2 1]);
        if (abs(xx(1)-x)>0.00001 || abs(xx(2)-n)>0.000001)
            keyboard
        end
        for k=1:i
            t=fscanf(fid2,'%f %f',[2 1]);
            Ryspolytable(j,i,k)=t(1);       %t
            Ryspolytable(j,i,k+i)=t(2);     %w
        end
    end
    end

    S=zeros(n_basis,n_basis);
    T=zeros(n_basis,n_basis);
    Va=zeros(n_basis,n_basis);
    Vb=zeros(n_basis,n_basis);
    V4=zeros(n_basis,n_basis,n_basis,n_basis);
    
    %calculate one-electron integrals
    %Overlap integral
    for k=1:n_basis  % the kth basis
    for l=k:n_basis  % the lth basis
        for i=1:n_pg_basis(k)
        for j=1:n_pg_basis(l)
            S(k,l)=S(k,l)+overlap_int(center_basis(k,:),alpha_basis(k,i),type_basis(k,:),center_basis(l,:),alpha_basis(l,j),type_basis(l,:))*coef_basis(k,i)*coef_basis(l,j);
            T(k,l)=T(k,l)+kinetic_energy_int(center_basis(k,:),alpha_basis(k,i),type_basis(k,:),center_basis(l,:),alpha_basis(l,j),type_basis(l,:))*coef_basis(k,i)*coef_basis(l,j);
            Va(k,l)=Va(k,l)-Za*NAI([0 0 0],center_basis(k,:),alpha_basis(k,i),type_basis(k,:),center_basis(l,:),alpha_basis(l,j),type_basis(l,:))*coef_basis(k,i)*coef_basis(l,j);
        end
        end
        S(l,k)=S(k,l);
        T(l,k)=T(k,l);
        Va(l,k)=Va(k,l);
    end
    end
    
    %calculate two-electron integrals
    for i=1:n_basis
    for j=i:n_basis
    for k=1:n_basis
    for l=k:n_basis
        for ii=1:n_pg_basis(i)
        for jj=1:n_pg_basis(j)
        for kk=1:n_pg_basis(k)
        for ll=1:n_pg_basis(l)
            V4(i,j,k,l)=V4(i,j,k,l)+TWOE(center_basis(i,:),alpha_basis(i,ii),type_basis(i,:),center_basis(j,:),alpha_basis(j,jj),type_basis(j,:),center_basis(k,:),alpha_basis(k,kk),type_basis(k,:),center_basis(l,:),alpha_basis(l,ll),type_basis(l,:),Ryspolytable)*coef_basis(i,ii)*coef_basis(j,jj)*coef_basis(k,kk)*coef_basis(l,ll);
        end
        end
        end
        end
        V4(j,i,k,l)=V4(i,j,k,l);
        V4(i,j,l,k)=V4(i,j,k,l);
        V4(j,i,l,k)=V4(i,j,k,l);
    end
    end
    end
    end
    
    output(S,T,Va,V4,fid);
end

    function [twoe]=TWOE(center1,alpha1,type1,center2,alpha2,type2,center3,alpha3,type3,center4,alpha4,type4,Ryspolytable)
        L=sum(type1)+sum(type2)+sum(type3)+sum(type4);
        n=floor(L/2)+1;
        RA=(alpha1*center1+alpha2*center2)/(alpha1+alpha2);
        RB=(alpha3*center3+alpha4*center4)/(alpha3+alpha4);
        A=alpha1+alpha2;
        B=alpha3+alpha4;
        rou=A*B/(A+B);
        AB=RA-RB;
        D=rou*AB.*AB;
        t=zeros(n,1);
        W=zeros(n,1);
        j=floor(sum(D)/0.1)+1;
        if (j<=1001)
            for i=1:n
                t(i)=Ryspolytable(j,n,i)*(j-sum(D)/0.1)+Ryspolytable(j+1,n,i)*(sum(D)/0.1-j+1);
                W(i)=Ryspolytable(j,n,i+n)*(j-sum(D)/0.1)+Ryspolytable(j+1,n,i+n)*(sum(D)/0.1-j+1);
            end
        else
            [t,W]=RysPoly(sum(D),n,crit);
        end
        twoe=0;
        for i=1:n
            u=rou*t(i)^2/(1-t(i)^2);  %herer u is actually u^2;
            B00=1/2/(A+B)*t(i)^2;
            B10=1/2/A-B/2/A/(A+B)*t(i)^2;
            B01=1/2/B-A/2/B/(A+B)*t(i)^2;
            C00=(RA-center1)+B*(RB-RA)/(A+B)*t(i)^2;
            C00_prime=(RB-center3)+A*(RA-RB)/(A+B)*t(i)^2;
            G=exp(-D*t(i)^2)*sqrt(1-t(i)^2)*pi/sqrt(A*B); %G00
            Ix=exp(D(1)*t(i)^2)/sqrt(1-t(i)^2)*I(B00,B10,B01,C00(1),C00_prime(1),G(1),type1(1),type2(1),type3(1),type4(1),center1(1),center2(1),center3(1),center4(1));
            Iy=exp(D(2)*t(i)^2)/sqrt(1-t(i)^2)*I(B00,B10,B01,C00(2),C00_prime(2),G(2),type1(2),type2(2),type3(2),type4(2),center1(2),center2(2),center3(2),center4(2));
            Iz=exp(D(3)*t(i)^2)/sqrt(1-t(i)^2)*I(B00,B10,B01,C00(3),C00_prime(3),G(3),type1(3),type2(3),type3(3),type4(3),center1(3),center2(3),center3(3),center4(3));
            twoe=twoe+2*sqrt(rou/pi)*W(i)*Ix*Iy*Iz;
        end
    end
    function [I]=I(B00,B10,B01,C00,C00_prime,G11,n1,n2,n3,n4,x1,x2,x3,x4)
%         ngrid=1;
%         cutoff=sqrt(10/max([alpha1,alpha2]));
%         x1min=-cutoff;
%         x1max=cutoff;
%         cutoff=sqrt(10/max([alpha3,alpha4]));
%         x2min=-cutoff;
%         x2max=cutoff;
%         delI=1.0;
%         I=-1;
%         while (delI>crit)
%             if (ngrid<100)
%                 ngrid=ngrid*2+1;
%             else
%                 ngrid=ngrid+100;
%             end
%             I0=I;
%             I=cal_I(n1,n2,n3,n4,x1,x2,x3,x4,alpha1,alpha2,alpha3,alpha4,u,x1min,x1max,x2min,x2max,ngrid);
%             delI=abs(I-I0);
%         end
%         h1=(x1max-x1min)/ngrid;
%         h2=(x2max-x2min)/ngrid;
%         delI=1.0;
%         while (delI>crit)
%             ngrid=floor(1.1*ngrid/2)*2+1;
%             x1min=-h1*ngrid/2;
%             x1max=h1*ngrid/2;
%             x2min=-h2*ngrid/2;
%             x2max=h2*ngrid/2;
%             I0=I;
%             I=cal_I(n1,n2,n3,n4,x1,x2,x3,x4,alpha1,alpha2,alpha3,alpha4,u,x1min,x1max,x2min,x2max,ngrid);
%             delI=abs(I-I0);
%         end
        %%% Another method, see J. Comp. Chem., 4, 154, 1983.
        %first, calculate all Gij. i=0,...,ni+nj; j=0,...,nk+nl.
        G=zeros(n1+n2+1,n3+n4+1);
        G(1,1)=G11;
        G(2,1)=C00*G(1,1);                           %G10
        for i=2:n1+n2                                %G(i,0)
            G(i+1,1)=(i-1)*B10*G(i-1,1)+C00*G(i,1);  
        end
        G(1,2)=C00_prime*G(1,1);                     %G01
        for i=2:n3+n4                                %G(0,j)
            G(1,i+1)=(i-1)*B01*G(1,i-1)+C00_prime*G(1,i);
        end
        for i=1:n1+n2
            G(i+1,2)=i*B00*G(i,1)+C00_prime*G(i+1,1);
        for j=2:n3+n4
            G(i+1,j+1)=(j-1)*B01*G(i+1,j-1)+i*B00*G(i,j)+C00_prime*G(i+1,j);
        end
        end
        %second, calculate Ix
        I2=0.0;
        for i=0:n2
        for j=0:n4
            I2=I2+c(i,n2)*(x1-x2)^i*c(i,n4)*(x3-x4)^j*G(n1+n2-i+1,n3+n4-j+1);
        end
        end
        I=I2;
%         if (((I/I2)>1.01 || (I/I2)<0.99) && abs(I)>0.000001)
%         [alpha1,alpha2,alpha3,alpha4,ngrid,I,I2]
%         end
    end
    function [I]=cal_I(n1,n2,n3,n4,x1,x2,x3,x4,alpha1,alpha2,alpha3,alpha4,u,x1min,x1max,x2min,x2max,ngrid)
        h1=(x1max-x1min)/ngrid;
        xv1=linspace(x1min+h1/2,x1max-h1/2,ngrid);
        h2=(x2max-x2min)/ngrid;
        xv2=linspace(x2min+h2/2,x2max-h2/2,ngrid);
        I=0;
        for ii=1:ngrid
        for jj=1:ngrid
            Q=alpha1*(xv1(ii)-x1)^2+alpha2*(xv1(ii)-x2)^2+alpha3*(xv2(jj)-x3)^2+alpha4*(xv2(jj)-x4)^2+u*(xv1(ii)-xv2(jj))^2;
            I=I+(xv1(ii)-x1)^n1*(xv1(ii)-x2)^n2*(xv2(jj)-x3)^n3*(xv2(jj)-x4)^n4*exp(-Q);
        end
        end
        I=I*h1*h2;
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
        j2=j;
        if (j2>l1)
            j2=l1;
        end
        for i=j1:j2    %in the term x^j, j=i+(j-i), i is from the 1st, (j-i) is from the 2nd
            f=f+c(i,l1)*a^(l1-i)*c(j-i,l2)*b^(l2-j+i);
        end
    end
    function [c]=c(i,j)  %calculate C^i_j
        c=1;
        for l=1:i
            c=c*(j-l+1)/l;
        end            
    end
    function [t]=kinetic_energy_int(center1,alpha1,type1,center2,alpha2,type2)
        l2=type2(1);
        m2=type2(2);
        n2=type2(3);
        t = alpha2*(2*(l2+m2+n2)+3)*overlap_int(center1,alpha1,type1,center2,alpha2,type2);
        type2temp=type2+[2 0 0];
        t = t - 2*alpha2^2*overlap_int(center1,alpha1,type1,center2,alpha2,type2temp);
        type2temp=type2+[0 2 0];
        t = t - 2*alpha2^2*overlap_int(center1,alpha1,type1,center2,alpha2,type2temp);
        type2temp=type2+[0 0 2];
        t = t - 2*alpha2^2*overlap_int(center1,alpha1,type1,center2,alpha2,type2temp);
        if (l2>1)
            type2temp=type2+[-2 0 0];
            t = t - 0.5*l2*(l2-1)*overlap_int(center1,alpha1,type1,center2,alpha2,type2temp);
        end
        if (m2>1)
            type2temp=type2+[0 -2 0];
            t = t - 0.5*m2*(m2-1)*overlap_int(center1,alpha1,type1,center2,alpha2,type2temp);
        end
        if (n2>1)
            type2temp=type2+[0 0 -2];
            t = t - 0.5*n2*(n2-1)*overlap_int(center1,alpha1,type1,center2,alpha2,type2temp);
        end
    end
    function [va]=NAI(centerN,center1,alpha1,type1,center2,alpha2,type2)
        l1=type1(1);
        m1=type1(2);
        n1=type1(3);
        l2=type2(1);
        m2=type2(2);
        n2=type2(3);
        gamma=alpha1+alpha2;
        P_C=(alpha1*center1+alpha2*center2)/(alpha1+alpha2);
        PA=P_C-center1;
        PB=P_C-center2;
        P=P_C-centerN;
        AB=center1-center2;
        AB2=AB*AB';
        PC2=P*P';
        imax=l1+l2;
        jmax=m1+m2;
        kmax=n1+n2;
        W=0;
        for mu=0:imax+jmax+kmax
            coef_F=0;
            for i=max(mu-jmax-kmax,0):min(mu,imax)
                for j=max(mu-i-kmax,0):min(mu-i,jmax)
                    k=mu-i-j;
                    coef_F=coef_F+G(i,l1,l2,PA(1),PB(1),P(1),gamma)*G(j,m1,m2,PA(2),PB(2),P(2),gamma)*G(k,n1,n2,PA(3),PB(3),P(3),gamma);
                end
            end
            W=W+coef_F*aug_int(PC2*gamma,mu);
        end
        va = 2*pi/gamma*exp(-alpha1*alpha2*AB2/gamma)*W;
    end
    function [g]=G(i,l1,l2,a,b,c,gamma)
        g=0;
        for ii=0:l1+l2
        for rr=0:floor(ii/2)
        for uu=0:floor((ii-2*rr)/2)
            if (ii-2*rr-uu == i)
                g=g+(-1)^(ii+uu)*f(ii,l1,l2,a,b)*factorial(ii)*c^(ii-2*rr-2*uu)*(1/4/gamma)^(rr+uu)/factorial(rr)/factorial(uu)/factorial(ii-2*rr-2*uu);
            end
        end
        end
        end
    end
    function [F]=aug_int(x,n)
    % This program computes the integral
    % int[t^2n * exp(-x*t^2)]
    % to the accuracy of crit
        crit=1e-5;
        delI=100;
        ngrid=1;
        a=0;
        b=1;
        F1=100;
        while (abs(delI-1) > crit)
            F0=F1;
            ngrid=ngrid*2+1;
            F1=intFgrid(x,n,a,b,ngrid);
            delI=F1/F0;
        end
        F=F1;
        
    end    
    function [F]=intFgrid(x,n,a,b,ngrid)

        h=(b-a)/ngrid;
        t_vec=linspace(a+h/2,b-h/2,ngrid);
        F=0.0;
        for i=1:ngrid
            f0=t_vec(i)^(2*n)*exp(-x*t_vec(i)^2);
            f2=(2*n*(2*n-1)*t_vec(i)^(2*n-2)-8*n*x*t_vec(i)^(2*n)+t_vec(i)^(2*n)*(4*x^2*t_vec(i)^2-2*x))*exp(-x*t_vec(i)^2);
            F=F+f0*h+f2*h^3/6;
        end
    end
    
    function output(S,T,Va,V4,fid)
        %write S, T, Va, Vb, V4
        %center A is first atom. center B is second atom
        %origin is on center A
        %V12A = off-diagonal nuclear attraction to center a, etc.
        fprintf(fid,'S\n');
        temp=size(S);
        nsize=temp(1);
            for i=1:nsize
            for j=1:nsize
                fprintf(fid,'%10.6f',S(i,j));
            end
            fprintf(fid,'\n');
            end
            fprintf(fid,'T\n');
            for i=1:nsize
            for j=1:nsize
                fprintf(fid,'%10.6f',T(i,j));
            end
            fprintf(fid,'\n');
            end
            fprintf(fid,'Va\n');
            for i=1:nsize
            for j=1:nsize
                fprintf(fid,'%10.6f',Va(i,j));
            end
            fprintf(fid,'\n');
            end
%             fprintf(fid,'V4\n');
%             for i=1:nsize
%             for j=1:nsize
%             for k=1:nsize
%             for l=1:nsize
%                 fprintf(fid,'%10.6f',V4(i,j,k,l));
%             end
%             fprintf(fid,'\n');
%             end
%             end
%             end
    end