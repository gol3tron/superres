%Script to Reprocess the particle tracking Data

%file paths for input (res files) and output files

basepath=   'F:\Thermal SPT\Thermal SPT organized\POPC\POPC SPT cold redo\neg9C 10ms 205x200\tracking\';
savepath=   'F:\Thermal SPT\Thermal SPT organized\POPC\POPC SPT cold redo\neg9C 10ms 205x200\analysis\';
imagepath=  'F:\Thermal SPT\Thermal SPT organized\POPC\POPC SPT cold redo\neg9C 10ms 205x200\';
filetype=1;  %1=tif files 2=matlab files
%base parmeters
first=1; %first image series
last=205; %lastimage series
startframe=1;
endframe=200;

%Camera Calibrations:
cal=0.2255975;  %0.18765766; %0.2255975; %5-11-12 60x = 0.33370  %5-10-12 100X=0.18765766; %6-30-2011 calibration = 0.18765766;  %6-22-2011 calibration emccd microscope = 0.189056154  %old calibration emccd microscope = 0.14805;  %microns per pixel of the camera

%tracking parameters
goodenough=9; %minimum number of frames from a track to be counted as a trajectory
memory=0; %number of frames a particle can drop out (blink)
tracknumber=''; %tracknumber
secperframe=0.01; %time associated with the individual frames in seconds
repeats=0; %minimum number of repeats at a particular stop
calculateontimes=0; %calculate on times? 1=yes 0=no
eliminatefirstframetracks=0; %eliminates tracks that begin on the first frame of the movie yes=1 no=0
calculateofftimes=0; %calculate off times? 1=yes 0=no
simpledecay=0; %calculate a simple decay? 1=yes 0=no
firstframeonly=0; %track particles that appear in first frame only 1=yes 0=no
minframefilter=0; %track particles that appear sometime after the first frame yes>0, no = 0
stdeviationfilter=1; %use a standard deviation from the mean filter to eliminate immobile particals 0=no, 1=x and y std, 2=r std
lowstd=0;  %200.0; %the low cutoff used in the standarddeviation filter (x, y corrdinates)
highstd=10000; %the high cutoff used in the standarddeviation filter (x, y corrdinates)
iMSDfilteroption=1;% was 2 9-23 %0=no filter, 1=max average MSD filter, 2= average MSD filter at a time lag
maxiMSD=0;%0.001 9-22-2011  %0.075
lagtime=6; %lagtime for with filtering takes place
iMSD_tail_slope=0.0; %0.015%set min slope of tail
iMSD_tail_intercept=0.00;
corralsizesquared=0.00; %filters by corralsizesquared
maxcs=10000000;
maxiSD=0.000;

%Calculate individual MSD and individual square dusplacements
calISD=0; %calcualte the individual squared displacements of each tracked particle 1=yes 0=no
calIMSD=1; %calculate the individual mean square displacements of each tracked particle
maxlag=6; %has to be less than or equal to isgood
rsq=0.5; %0.5; %sets rsquare value for goodness of fit, needs to be between 1 and 0; .90 says 90% of data can be explained by a linear model

%Diffusion Coefficient Parameter calculated from iMSD data
calcdiffusioncoefficent=1; % calculate diffision coefficient?  1=yes 0=no
GOF=0;%0.98 is good --0.91 %automatically chooses the model type depending on the rsquared value; yes>0
manualyevaluatetracktype=0; % manualy identify the type of diffusion for each track 1=yes 0=no
modeltype=1; %choose the model you wish to analyze data with 1=normal 2=anomalous 3=corralled
dplot=0; %plot difussion fits
keepall=1; %Choose what to keep=0 Keep all the fits =1
trackplot=0; %plot particle track?  1=yes, 0=no
plottime=1; %length of time you what fit if iMSD's vs. timelag to be displayed
bins=10;
bin_resolution=0.1;

%Displacement from mean
caldispcorrd=0; %calculate displacement corrdinates? 1=yes 0=no
DFMedge1=-400;
DFMedge2=400;
DFMresolution=0.0005; %in microns
type=2; %1=displacement from initial point 2=displacement from the mean

%Frame-to-Frame Step Size Distribution Analysis
SSD=0; %Do the step-size distribution analysis? 1=yes, 0=no
SSDedge1=0;
SSDedge2=200;
SSDresolution=0.002; %in microns

%animation paramerters
pl=0; %make an animation of tracks 1=yes 0=no
startplot=startframe; %designates start range of animation
endplot=endframe; %disginate the end range for the animation
CLOW=75; %color map for animation, low contrast setting
CHIGH=230; %color map for animation, high constrast setting
bandpass=1; %use the bandpass filter in plot? 1 for yes or 0 for no

%user eliminated tracks
track1=0;
track2=0;
track3=0;
track4=0;
track5=0;

%Finite Impulse Finction for Smoothing Data
lamda=1; %0.1; %average length scale of noise in pixels
w=3; %A length in pixels somewhat larger than *half* a typical object. Must be an odd valued integer.


cropmovie=0;
XMIN=55;
XMAX=70;
YMIN=30;
YMAX=45;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=first:last
    
%%%%%%%%%%track / diffusion constant program%%%%%%%%%%%%%%%%
    
    %load([basepath 'tracking\res_files\res_movie' num2str(j) '.mat']); %old
    load([basepath 'res_movie' num2str(j) '.mat']); %new

    if size(res,1)>0
        res=blink(res, memory);
    end
    
    if size(res,1)>0
        res=goodenoughframes_c(res, goodenough); %for a 9-coulmn results file
    end
    
    %if size(res,1)>0
        %res=[res res(:,8)];
    %end

    
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
            %stddev=standarddeviation(res,cal); %for x y std filtering
            %res=standarddeviationfilter_b(res,stddev,lowstd,highstd); %for
            %x y filtering
            %stddev=standarddeviation_b(res,cal); %for r filtering
            %res=standarddeviationfilter_radius_b(res,stddev,lowstd,highstd); %for r filtering
            %stddev=maximumdisplacement(res,cal); %for frame to frame r diplacement filtering
            %res=displacementfilter_with_range_radius(res,stddev,lowstd,highstd); %for r displacement filtering
            %stddev=standarddeviation_c(res,cal); %for x y r frame to frame std filtering
            %res=standarddeviationfilter_radius_b(res,stddev,lowstd,highstd); %for r filtering
            stddev=standarddeviation_from_mean_position_b(res,cal,goodenough); %for x y r std from mean position filtering
            res=standarddeviationfilter_with_options(stdeviationfilter,res,stddev,lowstd, highstd); %for either x and y, or r filtering
        end
    end
    
    if trackplot==1  %this both plots the track and sorts the track by type of diffusion test=lowstd
    if res>0 
        res=plot_track_plus_identify_results_e(res,XMIN,XMAX,YMIN,YMAX,cal,lowstd,goodenough,secperframe,rsq,maxlag,plottime,manualyevaluatetracktype,modeltype);
    end
    end
    
    %save( [savepath 'res_movie' num2str(j) '.mat'], 'res' );
    %%%%%%%%%%%%%%%%end of filtering for individual Movies%%%%%%%%%%%%%%
    
    
    
    if caldispcorrd==1 %does not use column 8 of res to identify track as mobile or immobile
   
        if size(res,1)>0
            dcorr=displacementcorrdinates_b(res,cal,type);
            save( [savepath 'dcorr' num2str(j) '.mat'], 'dcorr' );
            disp=histograms(dcorr,DFMedge1,DFMedge2,DFMresolution,goodenough);
            save( [savepath 'disp' num2str(j) '.mat'], 'disp' );
            clear disp dcorr;
        else
            dcorr=[];
            save( [savepath 'dcorr' num2str(j) '.mat'], 'dcorr' );
            disp=[];
            save( [savepath 'disp' num2str(j) '.mat'], 'disp' );
            clear disp dcorr;
        end
    
    end

    if SSD==1  %does not use column 8 of res to identify track as mobile or immobile
   
        if size(res,1)>0
            stepsize=frametoframestepsize(res,cal);
            save( [savepath 'stepsize' num2str(j) '.mat'], 'stepsize' );
            ssd=histograms(stepsize,SSDedge1,SSDedge2,SSDresolution,goodenough);
            save( [savepath 'ssd' num2str(j) '.mat'], 'ssd' );
            clear stepsize ssd;
        else
            stepsize=[];
            save( [savepath 'stepsize' num2str(j) '.mat'], 'stepsize' );
            ssd=[];
            save( [savepath 'ssd' num2str(j) '.mat'], 'ssd' );
            clear stepsize ssd;
        end
    
    end
    
    if calISD==1  %does not use column 8 of res to identify track as mobile or immobile
   
        if size(res,1)>0
            isd=individualSD_d(res,cal,secperframe,rsq,maxiSD,maxlag,dplot);
            save( [savepath 'isd' num2str(j) '.mat'], 'isd' );
            clear isd;
        else
            isd=[];
            save( [savepath 'isd' num2str(j) '.mat'], 'isd' );
            clear isd;
        end
    
    end
    
    if calIMSD==1  %uses column 8  of res (from standardeviation filtering) to identify track as mobile or immobile
   
        if size(res,1)>0
            iMSD=individualMSD_h(res,cal,secperframe);
            
            if iMSDfilteroption>0
                iMSD=iMSDfilter(iMSD,maxiMSD,iMSDfilteroption,lagtime);
                
            end
            
            %iMSD=individualMSD_c(res,cal,secperframe);
            save( [savepath 'iMSD' num2str(j) '.mat'], 'iMSD' );
            %clear iMSD;
        else
            iMSD=[];
            save( [savepath 'iMSD' num2str(j) '.mat'], 'iMSD' );
            %clear iMSD;
        end
        
        if size(iMSD,1)>0
            res=iMSDfilter_for_results_file(res,iMSD,goodenough,cal);   %res=iMSDfilter_for_results_file_b(res,iMSD,goodenough,cal) is for IGOR ploting with origin equal to zero and units in microns
        else
            res=[];
        end
        %save( [savepath 'res_movie' num2str(j) '.mat'], 'res' )
    end
    
    if calcdiffusioncoefficent==1  %diffusion_coefficent_linear_model used the iMSD files
        if size(iMSD,1)>0
            
            if GOF==0
                if manualyevaluatetracktype==1
                    diffusion=diffusion_coefficent_linear_anomalous_corralled_models(iMSD,maxlag,rsq,dplot,plottime);
                else
                    diffusion=diffusion_coefficent_choose_linear_anomalous_corralled(iMSD,maxlag,rsq,dplot,plottime,modeltype);
                    iMSD=reconcile_iMSD_file(iMSD,diffusion);
                    res=reconcile_res_file(res,diffusion);
                    if size(res,1)>0
                        quasi_D=quasi_diffusion_coefficients(res,cal,bins,bin_resolution,secperframe,dplot);
                        diffusion=[diffusion quasi_D];
                    end
                        save( [savepath 'iMSD' num2str(j) '.mat'], 'iMSD' );
                end
            else
                diffusion=diffusion_coefficent_linear_anomalous_corralled_e(iMSD,maxlag,dplot,GOF,corralsizesquared,iMSD_tail_slope,maxcs,maxiMSD,keepall);
                %diffusion=diffusion_coefficent_linear_or_corralled_b(iMSD,maxlag,dplot,GOF,corralsizesquared,iMSD_tail_slope,iMSD_tail_intercept); %9-23 0.015 (iMSD,maxlag,dplot,GOF,0.02)%note: diffusion=[D corral_size GOF_rsquare particle_type particle_identification_number];
                iMSD=reconcile_iMSD_file(iMSD,diffusion);
                res=reconcile_res_file(res,diffusion);
                %quasi_D=quasi_diffusion_coefficients(res,cal,bins,bin_resolution,secperframe);
                if size(res,1)>0
                quasi_D=quasi_diffusion_coefficients(res,cal,bins,bin_resolution,secperframe,dplot);
                diffusion=[diffusion quasi_D];
                end
                save( [savepath 'iMSD' num2str(j) '.mat'], 'iMSD' );
            end
            
            save( [savepath 'diffusion' num2str(j) '.mat'], 'diffusion' );
            clear diffusion;
            %old fitting rutine %trsort=[res(:,1:2) res(:,6) res(:,8)];
            %old fitting rutine %diffusion=diffco_MSD_slope_v4(trsort,secperframe,maxlag,rsq,dplot);
            %old fitting rutine %save( [savepath 'diffusion' num2str(j) '.mat'], 'diffusion' );
            %old fitting rutine %clear diffusion;
        else
            diffusion=[];
            save( [savepath 'diffusion' num2str(j) '.mat'], 'diffusion' );
            clear diffusion;
        end
    end
    clear iMSD;
    

    
%%%%%%%%%%%%%%%%%%%%%%%Animation Routine%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pl==1
    
   %%%% sets figure parameters %%%%
    CLIM=[CLOW CHIGH];
    colormap(gray);
    if bandpass==0
        if filetype==1
            strnam=[imagepath,'0',num2str(j),'.tif'];
            img=imread(strnam,1);
        end
        if filetype==2
            strnam=[imagepath,'crop',num2str(j),'.mat'];
            load([imagepath, 'crop' num2str(j) '.mat']);
            img=sim(:,:,1);
        end
    end
        
     if bandpass==1
        if filetype==1
            strnam=[imagepath,'0',num2str(j),'.tif'];
            image_data=imread(strnam,1);
        end
        if filetype==2
            strnam=[imagepath,'crop',num2str(j),'.mat'];
            load([imagepath, 'crop' num2str(j) '.mat']);
            image_data=crop(:,:,1);
        end
        img=bpass_brozik(image_data,lamda,w);
     end
    imagesc(img(:,:),CLIM);
    if cropmovie==1
    axis([XMIN XMAX YMIN YMAX]); %for Book
    else
    axis equal
    axis tight
    end
    hold all
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
     for i=startplot:endplot
        
        if bandpass==0
            if filetype==1
                strnam=[imagepath,'0',num2str(j),'.tif'];
                img=imread(strnam,i);
            end
            if filetype==2
                strnam=[imagepath,'crop',num2str(j),'.mat'];
                load([imagepath, 'crop' num2str(j) '.mat']);
                img=crop(:,:,i);
            end
        end
        
        if bandpass==1
            if filetype==1
                strnam=[imagepath,'0',num2str(j),'.tif'];
                image_data=imread(strnam,i);
            end
            if filetype==2
                strnam=[imagepath,'crop',num2str(j),'.mat'];
                load([imagepath, 'crop' num2str(j) '.mat']);
                image_data=crop(:,:,i);
            end
            img=bpass_brozik(image_data,lamda,w);
        end
        
        
        imagesc(img(:,:),CLIM);
        
        if cropmovie==1
            axis([XMIN XMAX YMIN YMAX]); %for Book
        else
            axis equal
            axis tight
        end
        
        hold all;
        
        
        for row=1:size(res,1)
            if res(row,6)==i
                %plot(res(row,1),res(row,2),'yo-','LineWidth',1.1,'MarkerSize',100); %for book
                %plot(res(row,1),res(row,2),'b.-','LineWidth',1.1,'MarkerSize',7); %for book
                plot(res(row,1),res(row,2),'yo-','LineWidth',1.1,'MarkerSize',30);
                plot(res(row,1),res(row,2),'b.-','LineWidth',1.1,'MarkerSize',5); 
                text(res(row,1),res(row,2), ['\leftarrow ',int2str(res(row,9))],'Color','y','FontSize',18);
            end
        end
    
        pause(0.1)
        animation(i-startplot+1)=getframe(gca);
        hold off;
    
     end
     
    close all;
    eval(['animation', int2str(j), '=animation;']);
    save([savepath 'animation' num2str(j) '.mat'], 'animation');
    clear animation;
    eval(['clear animation', int2str(j)]);

end

if size(res,1)>0
    res=shift_tracks_to_origin(res,cal);
    
    if trackplot==1
        plot_track(res);
    end
    
end

save( [savepath 'res_movie' num2str(j) '.mat'], 'res' );

end

clear all