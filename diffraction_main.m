%% Input paramters
tic
%Define the size of the extended source and initialise the variables
s_size=0e-6;
N_r=50;
N_theta_1=7;
res_r=s_size/N_r;

%Defining constants for the problem
lam=5e-11; %wavelength
D=15e-2; %Source to plate distance (13.4e-2 from skimmer to pinhole)
a_s=1;%sqrt(4.39e19); %Source amplitude

%Zone plate dimensions
f=2.12e-3; %Focal Point distance
max_r=3e-6; %Max dimension of aperture
min_feature=50e-9; %Min feature size

%Define the area over which the intensity is calculated
L_min=1.2e-3;   %Minimum z distance from zone plate
L_max=3.2e-3;   %Minimum z distance from zone plate
N_L=25;         %Number of slices in z to calculate intensity on
N=1000;         %Number of pixels to define in x and y directions for problem, actual value will be the nearest power of 2 bigger than this number
max_b=4.8e-6;   %Length of grid intensity is defined on

%Find the size of the zone plate you have defined
n_max_r=zone_size(lam,f,max_r,min_feature);
max_size=n_max_r*1.01;

%Calculate how many source points you are going to use to esimate how long the
%code will take extremely crudely
% N_tot=N_theta_1*N_r*(N_r+1)/2;
% time=7; %Estimate of how long each step will take, needs adjusting manually
% tot_time=N_tot*time;

%% Generate variables

L=linspace(L_min,L_max,N_L);

if max_b<max_size
    warning(max_b<max_size)
end


%Calculate the spacing size between the pixels
delta=(2*(max_b))/(N-1);


xa=-max_size:delta:max_size;
xb=-max_b:delta:max_b;

%Use an overall matrix with a size that is a power of 2, makes the fft faster.
Na=size(xa);
Nb=size(xb);
padded_size=2.^nextpow2(Nb(2));
clipped_num=floor((padded_size/2+Na(2)-Nb(2)/2)); %Only centre of the output matrix is used since convolution is only valid there
pad_size=padded_size-(2*clipped_num+1);

%Find x and y vectors
x_out=(delta*((1-pad_size/2):pad_size/2))';

%Generate matrix for the intensity
tot=zeros(pad_size,pad_size,N_L);
counter=1;
width=zeros(N_r,1);

%% Calculate aperture function and b matrix for problem

fftb=b_generation(L_min,L_max,N_L,N,max_size,max_b,lam);    %Generate b matrix and reuse. Uses more memory but saves on cpu time
h=aperture_funct(xa,xa,lam,n_max_r,f);  %Generate aperture function for zone plate
cut_surf=zeros(pad_size,N_L,N_r);
N_counter=0;

%% Main code
if s_size~=0    %Check the source is an extended source
    %Main loop over each radius for the source
    for r=res_r:res_r:s_size
        %Increasing number of points as you get further from centre of extended source to maintain roughly equal density
        N_theta=ceil(N_theta_1*r/res_r);
        theta=linspace(0,2*pi,N_theta);
        %x and y co-ordinates of the source points
        x0=r*sin(theta);
        y0=r*cos(theta);
        
        %Loop over each source point at this radius adding up the square of
        %the calculated amplitude
        for n=1:N_theta-1
            tot=tot+(abs(diffract_conv_new(L_min,L_max,N_L,N,max_size,max_b,lam,D,a_s,n_max_r,f,x0(n),y0(n),fftb,h))).^2;
            N_counter=N_counter+1;
        end
        
        %Find maximum of the intensity
        [maxval, maxloc] = max(tot(:));
        [maxloc_row, maxloc_col, maxloc_l] = ind2sub(size(tot), maxloc);
        
        %Take the slice with the maximum intensity
        temp2=reshape(tot(:,maxloc_col,:),[pad_size,N_L]);
        cut_surf(:,:,counter)=temp2;
        %         figure
        %         surf(L,x_out,abs(temp2))
        %         shading flat
        %         view(2)
        
        %Find width of the central spot
        transv=reshape(tot(:,maxloc_col,maxloc_l),[pad_size 1]);
        width(counter)=fwhm_fit(x_out,transv);
        
        %Increase counter for loop
        counter=counter+1;
        disp(['N = ' int2str(counter-1)])
        
        %Option to save all the variables at each loop
        %         vars=who;
        %         vars(strcmp(vars, 'fftb'), :) = [];
        %         save(['/tmp/Diffraction_Calculations/output' int2str(counter-1) '.mat'],vars{:});
        toc
    end
    
%If the source is just a point source the code below runs (for debugging normally)  
elseif s_size==0
    x0=0;
    y0=0;
    N_counter=1;
    %Calculate the intensity
    tot=(abs(diffract_conv_new(L_min,L_max,N_L,N,max_size,max_b,lam,D,a_s,n_max_r,f,x0,y0,fftb,h))).^2;
    
    %Calculate the maximum point
    [maxval, maxloc] = max(tot(:));
    [maxloc_row, maxloc_col, maxloc_l] = ind2sub(size(tot), maxloc);
    
    temp2=reshape(tot(:,maxloc_col,:),[pad_size,N_L]);
    cut_surf(:,:,counter)=temp2;
%     figure
%     surf(L,x_out,abs(temp2))
%     shading flat
%     view(2)
        
    transv=reshape(tot(:,maxloc_col,maxloc_l),[pad_size 1]);
    
    width(counter)=fwhm_fit(x_out,transv);
    toc
else
    disp('What have you done?')
end

%% Post processing of results
new_size=size(tot);
r=res_r:res_r:s_size;

%Normalising
open_area=sum(sum(aperture_funct(xa,xa,lam,n_max_r,f)))*delta^2;
tot_I=(a_s/D)^2*open_area;
tot=tot/N_counter;

%Spot Efficiency
transv_slice=tot(:,:,maxloc_l);
max_eff_r=max_b-2*max_size;
N_eff=floor(max_eff_r/delta);
I_eff=zeros(N_eff,1);
mat_x=((1:new_size(1))-maxloc_row).^2;
mat_y=((1:new_size(2))-maxloc_col).^2;

R=sqrt(bsxfun(@plus,mat_y,mat_x'));

% for k=1:N_eff
%     for i=1:new_size(1)
%         for j=1:new_size(2)
%             if R(j,i)<k^2+1e-4
%                 I_eff(k)=I_eff(k)+tot(i,j,maxloc_l);
%             end
%         end
%     end
% end

%Find the intensity in a slowly increasing radius from the maximum
for k=1:N_eff
    locs=R <(k+1e-4);
    I_eff(k)=sum(sum(transv_slice(locs)));
end



eff_r=(1:N_eff)*delta;
I_eff=(I_eff*delta^2)/(tot_I);
figure;plot(eff_r,I_eff,'b')

toc