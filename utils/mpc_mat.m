function [Lx,Lu]=mpc_mat(A,B,N)
    n=size(A,1);
    m=size(B,2);
    Lx=zeros(n*N,n);
    Lu=zeros(n*N,m*N);
    for i=1:N
        Lx((i-1)*n+1:i*n,:)=A^i;
        for j=1:i
            Lu((i-1)*n+1:i*n,(j-1)*m+1:j*m)=A^(i-j)*B;
        end
    end
end

    
    
