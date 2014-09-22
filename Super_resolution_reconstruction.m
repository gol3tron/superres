% Script to take Super Resolution data and make a Reconstructed Image
% Brozik & Goler, last update 14 August 2014

% definte the factor that converts pixel distances to real distances in
% nanometers. this is measured via calibration of a stationary fluor on the
% camera. the units of the conversion factor are nm/pixel.
conversion=187.65766;

% This parameter defines the desired resolution of the calculation in
% nanometers, i.e. it is used to round coordinates of reconstructed image
resolution_of_reconstructed_image=1; %5 %10

% loads the variable "cat", which should represent concatenated particle
% tracking output (is this correct?)
input_data = load('cat.mat');
%cat=MT; % used for testing 08 13 2014

% The super_data variable represents superresolution feature coordinate
% data. Concatenated tracking data with units of pixels are converted to
% nanometers with the specified resolution.
super_data=round(input_data*conversion/resolution_of_reconstructed_image); 

% At the moment, the super-resolution calcuation required only the x,y
% coordinates of features, which are found in the first two columns of
% super_data. Keep only those two columns. The variable super_data_full
% represents the full variable.
super_data_old = super_data;
super_data=[super_data(:,1) super_data(:,2)];

% The command sortrows is used to sort the x,y coordinates. In order to
% speed up the counting loops, a region of interest (ROI) is chosen to limit the
% number of steps necessary to process the data. In order to isolate this
% window, x,y coordinates between a min and max are discriminated out of
% the variable super_data such that only the data within the region of
% interest remain. The variable super_data_arch represents the variable
% super_data prior to this x,y cutting.
% The ROI is defined by the following min_x,miny, and max_x,max_y
% coordinates.
super_data_arch=super_data;
min_x=4505;
max_x=6771;
min_y=13439;
max_y=15500;

% super_data is sorted WRT the x-coordinates in column 1. Doing so allows
% us to cut data outside the ROI.
super_data=sortrows(super_data,1);

% min_x and max_x define, respectively, the minimum and maximum values of
% the x coordinate that defines the ROI. The find command determines the
% indices corresponding to the desired coordinates, and the following
% assignment command of super_data selects the desired window in the x
% direction.
index_min_x=find(super_data(:,1)==min_x);
index_max_x=find(super_data(:,1)==max_x);
super_data = super_data(index_min_x(end):index_max_x(1),:);

% The same processess that is used to select the x coordinates of the ROI
% is used for the y coordinates.
super_data=sortrows(super_data,2);
index_min_y=find(super_data(:,2)==min_y);
index_max_y=find(super_data(:,2)==max_y);
super_data = super_data(index_min_y(end):index_max_y(1),:);

% The variable image represents a blank matrix that will contain occurances
% of fluorescence detection at the corresponding x,y coordinates. Every
% time a coordinate is found within the ROI, the value of image at that
% coordinate is incremented.
image = zeros(max_x-min_x,max_y-min_y);
%image=zeros(max_x,max_y); %Old

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
    
    keepingtrack=keepingtrack-1;
end    
end


image(1:(min_x-1),:)=[];
image(:,1:(min_y-1))=[];

maximagevalue=max(max(image(:,:,1)));
minimagevalue=min(min(image(:,:,1)));

grayscaleimage(image, minimagevalue, maximagevalue);
