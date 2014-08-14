%Script to take Super Resolution data and make a Reconstructed Image

conversion=187.65766;  %conversion factor from pixels to nanometers in nanometers / pixel
resolution_of_reconstructed_image=1; %set reconstruction resolution in nm


load('cat.mat')
%cat=MT;

super_data=round(cat*conversion/resolution_of_reconstructed_image); %does conversion from pixels to chosen resolution
super_data=[super_data(:,1) super_data(:,2)];

super_data=sortrows(super_data,1);

super_data_arch=super_data;

min_x=find(super_data(:,1)==4505);
max_x=find(super_data(:,1)==6810);


super_data = super_data(min_x(end):max_x(1),:);

super_data=sortrows(super_data,2);

min_y=find(super_data(:,2)==13439);
max_y=find(super_data(:,2)==15500);

super_data = super_data(min_y(end):max_y(1),:);

min_x=4505;
max_x=6771;
min_y=13439;
max_y=15500;

image=zeros(max_x,max_y);

keepingtrack=(max_y-min_y+1)*(max_x-min_x+1);
for a=min_x:max_x
for c=min_y:max_y
    
    counter=0;
    for b=1:size(super_data,1)
        if super_data(b,1)==a
        if super_data(b,2)==c    
            counter=counter+1;
        else
            counter=counter;
        end
        else
            counter=counter;
        end
    end
    image(a,c)=counter;
    
    keepingtrack=keepingtrack-1
end    
end


image(1:(min_x-1),:)=[];
image(:,1:(min_y-1))=[];




maximagevalue=max(max(image(:,:,1)));
minimagevalue=min(min(image(:,:,1)));

grayscaleimage(image, minimagevalue, maximagevalue);
