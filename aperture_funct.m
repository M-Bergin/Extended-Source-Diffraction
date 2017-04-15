function h=aperture_funct(x,y,lam,n_max_r,f)
%Function to define the aperture function for a fresnel zone plate with
%supporting bars

% N=length(x);
% M=length(y);
%
% square_dim=2e-7;
%
% h=zeros(N,M);
%
% for i=1:N
%     for j=1:M
%         if (x(i)<square_dim && y(j)<square_dim && x(i)>-square_dim && y(j)>-square_dim)
%             h(i,j)=1;
%         end
%     end
% end

bar_width=50e-9/2;
N=length(x);
M=length(y);
h=zeros(N,M);
for i=1:N
    for j=1:M
        if sqrt(x(i).^2+y(j).^2)<n_max_r
            if floor(mod((x(i).^2+y(j).^2)/(lam*f),2))==1
                if x(i)<bar_width && x(i)>-bar_width
                    h(i,j)=0;
                elseif y(j)<bar_width && y(j)>-bar_width
                    h(i,j)=0;
                else
                    h(i,j)=1;
                end
            else
                h(i,j)=0;
            end
        else
            h(i,j)=0;
        end
    end
end

end
