function width=fwhm_fit(x,y)

maximum=max(y);

ind=find(y>maximum*0.2);

gaus_y=y(ind);
gaus_x=x(ind);

if(length(ind)<3)
    ind=find(y>maximum*0.15);

gaus_y=y(ind);
gaus_x=x(ind);
end

if(length(ind)<3)
    ind=find(y>maximum*0.1);

gaus_y=y(ind);
gaus_x=x(ind);
end

if(length(ind)<3)
    ind=find(y>maximum*0.05);

gaus_y=y(ind);
gaus_x=x(ind);
end

f=fit(gaus_x,gaus_y,'gauss1');
w=f.c1;

width=2*w;
end
