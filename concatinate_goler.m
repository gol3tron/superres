
savepath='/Volumes/Verbatim/Barden/Thermal Control c/DMPC/02072013/27C 10ms 200x205/analysis/';
file='diffusion';
startfile=1;
endfile=10;
columns=8;

cat=zeros(1,columns+1);
for i=startfile:endfile
     load([savepath file num2str(i) '.mat']);
     if size(diffusion)>0
     movienumber=i;
     mn=zeros(1,1);
     for ii=1:size(diffusion,1)
         mn=[mn;movienumber];
     end
     mn(1,:)=[];
     aa=[diffusion mn];
     cat=[cat; aa];
     end
end
cat(1,:)=[];

clear aa file i columns endfile startfile savepath diffusion res