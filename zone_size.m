function n_max_r=zone_size(lamd,f,max_r,min_feature)



n=1;

r(1)=0;

while r(n)<max_r
    r(n+1)=sqrt(n*lamd*f);
    n=n+1;
end
r(n)=[];

N=length(r);

for i=1:N-1
    if(r(i+1)-r(i)>min_feature)
        diff(i)=r(i+1)-r(i);
        cut=i;
    end
end

zone_rad=r(1:cut);
n_max_r=r(cut+1);

% figure;
% hold on;
% axis square;
% t = linspace(0,2*pi);
% for j=cut:-1:2
%     fill(zone_rad(j)*cos(t),zone_rad(j)*sin(t),[mod(j,2),mod(j,2),mod(j,2)]);
% end
% 
% 
% N_point=10000;
% x=linspace(-1.1*max_r,1.1*max_r,N_point);
% y=linspace(-1.1*max_r,1.1*max_r,N_point);
% h=zeros(N_point,N_point);
% 
% for i=1:N_point
%     for j=1:N_point
%         if sqrt(x(i).^2+y(j).^2)<n_max_r
%             if round(mod((x(i).^2+y(j).^2)/(lamd*f),2))==1
%                 h(i,j)=1;
%             end
%         end
%     end
% end

