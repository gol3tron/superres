%Script for particle tracking

%parameter file for the mpretrack_init and mpretrack_brozik programs
trtype=1; %track type initial setup=0; actual tracking=1;
basepath= 'D:\Thermal SPT\Thermal Control c\DMPC\01292013\27C 10ms 200x205\';
savepath= 'D:\Thermal SPT\Thermal Control c\DMPC\01292013\27C 10ms 200x205\tracking\';
filetype=1; %tiff data = 1 and matlab mov files = 2
prefix='crop'; % ='bin' for binned frames; ='sim' for simulation
bandpassfilter=1; % Do you want to use a bandpass and gaussian filter to track your data? 1=yes 0=use raw data in tracking program
binframes=0; %Note: if you bin frames you filetype=2 %Do you want to bin frames? 0=no and 0>yes with the number being the number of frame you wish to bin together

%base parmeters to pick points
featsize=3; %Radius of feature you are interested in; has to be an odd number
barint=1300;   %1300;  %minimum integrated intensity that will be accepted
barrg=4.5; %maximum radius of gyration squared that will be accepted
barcc=.3; %maximum eccentricity that will be accepted
IdivRg=10; %minimum ratio of integrated intensity to radius of gyration squared accepted
Imin=90;  %100; %minimum intensity of local maximum to be considered
masscut=100; %threshold for integrated intensity of features before position refinement
first=1; %first image series
last=10; %lastimage series
startframe=1; %201;
endframe=200; %300;
fovn=1; %numeridentifying the the series of images you are analyzing for intital setup
frame=100; %frame number in which you are optimizing parameters
field=2; %set to 0 or 1 for one field of an interlaced frame and 2 if it if a full frame

%Camera Calibrations:
cal=0.5; %0.18765766; %60x objective 5-11-12 = 0.33370  %5-10-12 100X=0.18765766; %6-30-2011 calibration = 0.18765766;  %6-22-2011 calibration emccd microscope = 0.189056154  %old calibration emccd microscope = 0.14805;  %microns per pixel of the camera

%tracking parameters
maxdisp=7; %7; %maxiumum displacement between frames in a track
goodenough=1; %2; %minimum number of frames from a track to be counted as a trajectory
memory=0; %number of frames a particle can drop out (blink)
secperframe=0.010; %time associated with the individual frames in seconds
eliminatefirstframetracks=0; %eliminates tracks that begin on the first frame of the movie yes=1 no=0
firstframeonly=0; %track particles that appear in first frame only 1=yes 0=no
minframefilter=0; %track particles that appear sometime after the first frame yes>0, no = 0
stdeviationfilter=0; %use a standard deviation from the mean filter to eliminate immobile particals 1=yes, 0=no
cutoffstd=100; %the cutoff used in the standarddeviation filter (x, y corrdinates)

%animation paramerters
pl=0; %make an animation of tracks 1=yes 0=no
startplot=startframe; %designates start range of animation
endplot=endframe; %disginate the end range for the animation
CLOW=50; %color map for animation, low contrast setting
CHIGH=200; %color map for animation, high constrast setting
bandpass=1; %use the bandpass filter in plot? 1 for yes or 0 for no

%user eliminated tracks
track1=0;
track2=0;
track3=0;
track4=0;
track5=0;

%user defined area of interest
areaofinterest=0; %area of interst 0=no, 1 = a single point to track, 2 = mutilple points to track 
singlespot=0; % are you tracking a single spot over many movies? 1=yes 0=no this corrects for any mechanical creep in the microscope
pixelrange=3; % range of pixels above and below center point
xpixel=67.594;
ypixel=17.812;

%Finite Impulse Finction for Smoothing Data
lamda=1;%1; %average length scale of noise in pixels
w=3;%3; %A length in pixels somewhat larger than *half* a typical object. Must be an odd valued integer.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
density=zeros(1,1);
%lengthoftrack=zeros(1,1);

xlow=xpixel-pixelrange;
xhigh=xpixel+pixelrange;
ylow=ypixel-pixelrange;
yhigh=ypixel+pixelrange;

if binframes>0
    Bin_and_Reconstruct_Movie_c(basepath,basepath,first,last,binframes);
    endframe=endframe/binframes;
    endplot=endframe;
    secperframe=secperframe*binframes;
    
    time=zeros(1,1);
    for bbb=1:endframe
        time=[time;bbb*secperframe];
    end
    time(1,:)=[];
    save( [basepath,'times.mat'], 'time' );
else
    time=zeros(1,1);
    for bbb=1:endframe
        time=[time;bbb*secperframe];
    end
    time(1,:)=[];
    save( [basepath,'times.mat'], 'time' );
end
clear time

if trtype==0
    %mpretrack_init_brozik_d(basepath, featsize, barint, barrg, barcc, IdivRg, fovn, frame, masscut, Imin, field, lamda, w, CLOW, CHIGH)
    mpretrack_init_brozik_bandpass_c(basepath, featsize, barint, barrg, barcc, IdivRg, fovn, frame, masscut, Imin, field, lamda, w, CLOW, CHIGH,filetype,prefix);
end


for j=first:last

%%%%%%%%%%track / diffusion constant program%%%%%%%%%%%%%%%%
if trtype==1

    mpretrack_brozik_g(basepath,savepath,j,barint,barrg,barcc,IdivRg,startframe,endframe,masscut,Imin,field, lamda, w, filetype,bandpassfilter,prefix);

    if areaofinterest>0
        if areaofinterest==1
            MT=spatiallimitset(j,savepath,xlow,xhigh,ylow,yhigh);
        end
        if areaofinterest==2
            MT=spatiallimitset_b(j,savepath,corrdinates,pixelrange);
        end
    end

    load([savepath 'MT_' num2str(j) '.mat']);

    if size(MT,1)>=goodenough
        if size(MT,1)>1
            fancytrack_brozik_d(basepath, savepath, j, maxdisp, 1, endframe-1);
        else
            res=[];
            save( [savepath 'res_movie' num2str(j) '.mat'], 'res' );
        end
    else
        res=[];
        save( [savepath 'res_movie' num2str(j) '.mat'], 'res' );
    end
    
    load([savepath 'res_movie' num2str(j) '.mat']);

    if size(res,1)>0
        res=blink(res, memory);
    end
    
    if size(res,1)>0
        res=goodenoughframes_b(res, goodenough); %for a 8-coulmn results file
    end
    
    if size(res,1)>0
    res=[res res(:,8)];
    end

    
    if size(res,1)>0
        if firstframeonly==1
            res=firstframefilter(res,1,1);
        end
    end
    
    if size(res,1)>0
        if minframefilter >0
            res=firstframefilter(res,0,minframefilter);
        end
    end
    
    if size(res,1)>0
        if stdeviationfilter>0
            stddev=standarddeviation_from_mean_position(res,cal,goodenough); %for x y r std from mean position filtering
            res=standarddeviationfilter_with_options(stdeviationfilter,res,stddev,lowstd, highstd); %for either x and y, or r filtering
        end
    end
    
    
    save( [savepath 'res_movie' num2str(j) '.mat'], 'res' );
    %%%%%%%%%%%%%%%%end of filtering for individual Movies%%%%%%%%%%%%%%
    
    
    

if pl==1
    
   %%%% sets figure parameters %%%%
    CLIM=[CLOW CHIGH];
    colormap(gray);
    if bandpass==0
        if filetype==1
            strnam=[basepath,'0',num2str(j),'.tif'];
            img=imread(strnam,1);
        end
        if filetype==2
            strnam=[basepath,prefix,num2str(fovn),'.mat'];
            load([basepath,prefix,num2str(fovn),'.mat']);
            eval(['img=',prefix,'(:,:,1);']);
        end
    end
        
     if bandpass==1
        if filetype==1
            strnam=[basepath,'0',num2str(j),'.tif'];
            image_data=imread(strnam,1);
        end
        if filetype==2
            strnam=[basepath,prefix,num2str(fovn),'.mat'];
            load([basepath,prefix,num2str(fovn),'.mat']);
            eval(['image_data=',prefix,'(:,:,1);']);
        end
        img=bpass_brozik(image_data,lamda,w);
     end
    imagesc(img(:,:),CLIM);
    axis equal;
    axis tight
    hold all
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
     for i=startplot:endplot
        
        if bandpass==0
            if filetype==1
                strnam=[basepath,'0',num2str(j),'.tif'];
                img=imread(strnam,i);
            end
            if filetype==2
                strnam=[basepath,prefix,num2str(fovn),'.mat'];
                load([basepath,prefix,num2str(fovn),'.mat']);
                eval(['img=',prefix,'(:,:,i);']);
            end
        end
        
        if bandpass==1
        if filetype==1
            strnam=[basepath,'0',num2str(j),'.tif'];
            image_data=imread(strnam,i);
        end
        if filetype==2
            strnam=[basepath,prefix,num2str(j),'.mat'];
            load([basepath,prefix,num2str(j),'.mat']);
            eval(['image_data=',prefix,'(:,:,i);']);
        end
        img=bpass_brozik(image_data,lamda,w);
        end
        
        
        imagesc(img(:,:),CLIM);
        axis equal;
        axis tight
        hold all;
        
        for row=1:size(res,1)
            if res(row,6)==i
                plot(res(row,1),res(row,2),'yo-','LineWidth',1.1,'MarkerSize',30);
                plot(res(row,1),res(row,2),'b.-','LineWidth',1.1,'MarkerSize',1); 
                text(res(row,1),res(row,2), ['\leftarrow ',int2str(res(row,8))],'Color','y','FontSize',18);
            end
        end
    
        pause(0.1)
        animation(i-startplot+1)=getframe(gca);
        hold off;
    
     end
     
    close figure 1;
    eval(['animation', int2str(j), '=animation;']);
    save([savepath 'animation' num2str(j) '.mat'], 'animation');
    clear animation;
    eval(['clear animation', int2str(j)]);

end


end



end


clear all
%clear calculateofftimes calculateontimes eliminatefirstframetracks simpledecay MT areaofinterest j pixelrange savepath singlespot xhigh xlow yhigh ylow xpixel ypixel off meanrepeatingspot track1 track2 track3 track4 track5 repeatingspot ii jj kk ll ontime ontimes spot sub subsub repeats on lngoftrack timelengthpertrack timestamp timestamps count particlenumber secperframe tracklength cal area den particlesperframe totframe tracknumber endframe startframe first last trtype lamda w pl IdivRg Imin barcc barint barrg basepath featsize field fovn frame goodenough masscut maxdisp memory numframes CHIGH CLOW CLIM i img row strnam bandpass image_data lownoise
