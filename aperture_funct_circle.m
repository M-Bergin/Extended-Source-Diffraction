function h=aperture_funct_circle(x,y,~,~,~)
%
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

circle_d=1.2e-6;
N=length(x);
M=length(y);
h=zeros(N,M);
for i=1:N
    for j=1:M
        if sqrt(x(i).^2+y(j).^2)<circle_d/2
            h(i,j)=1;
        else
            h(i,j)=0;
        end
    end
end
%figure;imagesc(h)
end
