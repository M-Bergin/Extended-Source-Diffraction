function fftb=b_generation(L_min,L_max,N_L,N,max_size,max_b,lam)
% Function to generate b matrix for diffraction calculation

%Find size of the spacing between points
delta=(2*(max_b))/(N-1);

%Generate vectors for matrix a
xa=-max_size:delta:max_size;
ya=-max_size:delta:max_size;

Na=size(xa);

%Generate vectors for matrix b
xb=-max_b:delta:max_b;
yb=-max_b:delta:max_b;

Nb=size(xb);

%Increase the size of the b matrix so that it is a power of 2
padded_size=2.^nextpow2(Nb(2));

fftb=zeros(padded_size,padded_size,N_L);

ik_vec= 1i*2*pi/lam;

L=linspace(L_min,L_max,N_L);

%Loop over each z position
for n=1:N_L
    
    %Calculate the b matrix
    b_t= L(n)* bsxfun(@(x2, y2)exp(ik_vec*(sqrt(L(n)^2+x2.^2+y2.^2)))./(L(n)^2+x2.^2+y2.^2), xb, yb');
    
    n_size_b=floor((padded_size-size(b_t))/2);
    
    %Pad the matrix with 0s to a power of 2
    b=padarray(b_t,n_size_b);
    
    if length(b)==padded_size-1
        b=[zeros(padded_size-1,1),b];
        b=[zeros(1,padded_size);b];
    end
    
    
    %IFFT
    fftb(:,:,n)=fft2(b);
    
end



