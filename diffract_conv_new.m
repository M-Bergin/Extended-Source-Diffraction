function result=diffract_conv_new(L_min,L_max,N_L,N,max_size,max_b,lam,D,a_s,n_max_r,f,x0,y0,fftb,h)
%Function to perform the convolution to calculate diffracted intensity from
%a defined aperture

% r2= @(x,y) (L^2+x.^2+y.^2); %actually this is r(x-xp,y-yp)
% s= @(x,y) sqrt(D^2+(x-x0).^2+(y-y0).^2);

% K= @(x,y) 0.5*(L./r(x,y));


%Pad to nearest power of 2 to increase fft efficiency
% pad_size=2.^nextpow2(N);


%Calcuate the vectors for the a and b matrices
delta=(2*(max_b))/(N-1);

xa=-max_size:delta:max_size;
ya=-max_size:delta:max_size;

Na=size(xa);

xb=-max_b:delta:max_b;
yb=-max_b:delta:max_b;

Nb=size(xb);

%Pad the b matrix to the next power of 2
padded_size=2.^nextpow2(Nb(2));

L=linspace(L_min,L_max,N_L);

%Calculate the a matrix
pre_fac=-(1i/lam)*delta^2*a_s;

a_t=pre_fac.*aperture_funct(xa,ya,lam,n_max_r,f).* ...
    bsxfun(@(x2, y2) exp(1i*2*pi.*(sqrt(D^2+(x2-x0).^2+(y2-y0).^2))/lam)./(sqrt(D^2+(x2-x0).^2+(y2-y0).^2)), xa, ya');

n_size_a=floor((padded_size-size(a_t))/2);


%Pad with zeroes the next power of 1
a=padarray(a_t,n_size_a);

if length(a)==padded_size-1
    a=[zeros(padded_size-1,1),a];
    a=[zeros(1,padded_size);a];
end

%Calculate the FFT
FFT1=fft2(a);

result=zeros(padded_size,padded_size,N_L);


for n=1:N_L
        
    %Correlate
    FFTR=FFT1.*fftb(:,:,n);
    
    %IFFT
    result(:,:,n)=fftshift(ifft2(FFTR)); %Can this be moved so it is only performed once? 
    
end
% figure;surf(x,y,abs(result))
% shading flat

%Find the middle of the result which is valid and return the result
%Rounding Errors Introduced here!!

clipped_num=floor((padded_size/2+Na(2)-Nb(2)/2));

result(end-clipped_num:end,:,:)=[];
result(1:clipped_num,:,:)=[];

result(:,end-clipped_num:end,:)=[];
result(:,1:clipped_num,:)=[];


