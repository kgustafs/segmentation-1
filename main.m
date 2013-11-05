%% experiment paramaters

cd /Users/bieler/Desktop/matlab/23september_40deg/

experimentPara;

% Define some stuff, move into the right folder

movie = 10;

outDir = [mainDir '/movie' num2str(movie) '/'];

cd(outDir)

N1=512;
N2=512;

addpath ../code

N = expe.Nframe;

%% combine stacks, denoise, ...

Nz = 3; %number of images to combine

deNoise = {'none','BM3D','median','localNorm'}; %denoise algo on each stack
deNoise = deNoise{3};

medianSize = 2;

weightsSegmentation = [0.8 0.9 1]; 
weightsData = [1 0.5 0.3];

doDraw = 0;

mkdirIfNotExist('zStackedYFP');
mkdirIfNotExist('zStackedYFP_Data');

for k=1:N

    disp(100*k/N);    
    
    prefix = [outDir 'img/' expe.filePrefix '_' num2str(k)  '_'];

    images = {};
    for i=1:Nz
        imagesPath{i} = [prefix num2str(i) '.tif'];
    end

    out = combineStack(imagesPath,N1,N2,Nz,deNoise,medianSize,weightsSegmentation,weightsData,doDraw);
    
    imwrite(out,['zStackedYFP/' num2str(k) '.png']);
    
    
    %just copy the raw images for quantification later
    
    out = imread(imagesPath{1});
    
    imwrite(out,['zStackedYFP_Data/' num2str(k) '.png']);
    
end

%% make a gif with the data

for k=1:N
    
        disp(100*k/N)
                 
        a = imread(['zStackedYFP/' num2str(k) '.png']);
        a = repmat(a,[1,1,3]);
        
        [a,cm] = rgb2ind(a,256);        
        if k == 1;
            imwrite(a,cm,['dataStacked.gif'],'gif', 'Loopcount',inf,'DelayTime',0.15);
        else
            imwrite(a,cm,['dataStacked.gif'],'gif','WriteMode','append','DelayTime',0.15);
        end
end

!open dataStacked.gif

%% just display the data

clf
colormap jet

m = [];

for k=1:N
    
        disp(100*k/N)
                 
        a = imread(['zStackedYFP/' num2str(k) '.png']);
        
        imagesc(a)
        caxis([0 150])
        pause(0.1)
        drawnow;
        
end

%% test segmentation on a few frames 
       
% inputFolder = 'zStackedYFP/';
% 
% testImages = [1 50  110 ];
% doDraw = 0;
% 
% clear tPara
% tPara(1).method = 'normal';
% 
% %! sizes are not interpolated in time
% 
% tPara(1).minSize= 1.7; %para of the first frame
% tPara(1).maxSize= 3.0;
% tPara(1).thresh = 0.91;
% tPara(1).openSize = 5.0;
% 
% tPara(2).minSize= 1.7; %para of the last frame
% tPara(2).maxSize= 3.0;
% tPara(2).thresh = 0.85;
% tPara(2).openSize = 5.0;
% 
% testSegmentationInterp;

edit guiSeg.m

%% Do the segmentation

inputFolder = 'zStackedYFP/';

load segPara.mat
load segParaFrames.mat

doDraw=0;

for frame=1:N
    disp(100*frame/N);
    
    para = [];
        
    para.minSize  = interp1(segParaFrames,[segPara(:).minSize  ],frame);
    para.maxSize  = interp1(segParaFrames,[segPara(:).maxSize  ],frame);
    para.thresh   = interp1(segParaFrames,[segPara(:).thresh   ],frame);
    para.openSize = interp1(segParaFrames,[segPara(:).openSize ],frame);
    para.ratioThresh = interp1(segParaFrames,[segPara(:).ratioThresh ],frame);
    
    %generateFilters;
        
    [filters openFilter] = generateFilters(para,doDraw);
    
    segmentImage(frame,para,inputFolder,filters,openFilter,0);
end

%% make segmentation gif 

threshFolder = 'zStackedThresh/';

makeSegmentationGif;
!open segmentation.gif

%% do the measures

threshFolder = 'zStackedThresh/';
doMeasures;

% try to split some merged cells

saveFolder = 'zStackedThreshSplit/';
mkdirIfNotExist(saveFolder)

doDraw = 0;

threshold = 1.5; %low value -> split everything

splitMergeCells;

threshFolder = saveFolder;
doMeasures;

%% look at difference between splited and original

doDraw = 1;
if doDraw 
    for k=1:N
        
         a = imread(['zStackedThresh/' num2str(k) '.png']);
         b = imread(['zStackedThreshSplit/' num2str(k) '.png']);

         a = double(a);
         b = double(b);

         imagesc(a+b);
         pause(0.1)
    end
end

%% correct segmentation by hand
% !first run the segmentation once

edit linksGui.m

%% redo the measures with corrected images

threshFolder = 'zStackedThreshCorrected/';
doMeasures;

%% load all measures & do Tracking

NToTrack = N;

Me = cell(1,N);

for k=1:N
    
    name = ['Measures/' num2str(k) '.mat'];
    load(name)
    
    tmp = [];
    for i=1:length(Measurements)
          tmp = [tmp Measurements(i)];
    end
    
    Me{k} = tmp;        
end

%define data given to the tracking algo
clear nuclei
for i=1:NToTrack    
   m=Me{i};
    
   nuclei(i).carth = cat(1,m.Centroid);
   nuclei(i).properties = [1/80*cat(1,m.MeanIntensity)];% 0+1*cat(1,m.Area)];
end

save nuclei.mat nuclei;

% do the tracking 

addpath ../code/Tracking/
addpath ../code/Tracking/libraries/

doDraw = 0; imageFraction = 1;

if( exist('zStackedThreshCorrected/1.png','file') )
    doLinksOnly = 0;
else    
    doLinksOnly = 1;
end

doTracking;

clf; imagesc(signal)

save tracks.mat tracks
save signal.mat signal
save traj.mat traj
save ind.mat ind
save divisions.mat divisions

%%

load tracks.mat 
load signal.mat 
load traj.mat
load ind.mat 
load divisions.mat 


%% plot Tracking

pauseTime = 0.05;

plotTracking;

%% select good traces, based on length

lengthOfTrace = zeros(size(ind,1),1);
lengthOfGaps  = zeros(size(ind,1),1);

indAnnotation = zeros(size(ind));

if( exist('lengthThresh.mat','file') )
    load lengthThresh.mat;
else
    lengthThresh = 0.5;
end

for i=1:size(ind,1)
        
    %look for continous traces
    indAnnotation(i,:) = markTrace(ind(i,:));

    lengthOfGaps(i) = sum( indAnnotation(i,:) == -3);
    lengthOfTrace(i) = sum( indAnnotation(i,:) == 1);
    
end

longTraces = find( (lengthOfGaps < 1) .* (lengthOfTrace/NToTrack>lengthThresh) );
clf

A = signal(longTraces,:);
[tmp,ia,ic] = unique(A,'rows');

longTraces = sort( longTraces(ia) );

imagesc(signal(longTraces,:))

length(longTraces)
colormap jet

save lengthThresh.mat lengthThresh
save longTraces.mat longTraces

%% build peak matrix

minTimeBetweenPeaks = 5;
peakMethod ='diff';
doDraw =0;

peakMatrix = zeros(size(ind)); 


for i=1:length(longTraces)
    
    idx = longTraces(i);
    tmp = signal(idx,:);

    maxp = getPeaks(tmp,expe,peakMethod,doDraw);
    
    for j=1:size(maxp,1)
        peakMatrix(idx,maxp(j,1))=1;
    end
    
end

%remove peaks that are too close in time    
for i=1:length(longTraces)
    
    idx = longTraces(i);
    tmp = signal(idx,:);
    
    pp = find(peakMatrix(idx,:));     
    
    for j=1:length(pp)-1

        p2pTime = expe.t(pp(j+1))-expe.t(pp(j));        
        if( p2pTime < minTimeBetweenPeaks)           

                %keep the best one
                score1 = tmp(pp(j));                
                score2 = tmp(pp(j+1));

                if(score1 >= score2)                   
                    peakMatrix(idx,pp(j+1)) = 0;
                else
                    peakMatrix(idx,pp(j)) = 0;
                end

        end
    end
end


%peaks at the begining are almost always false

peakMatrix(:,1:2) = 0;

imagesc(peakMatrix(longTraces,:) )
save peakMatrix.mat peakMatrix


%% build division matrix

minTimeBetweenDivisions = 4; 

divMatrix = zeros(size(ind)); 

for i=1:length(longTraces)
    
    idx = longTraces(i);
    tmp = signal(idx,:);
    
    %tracking division
    dd = find([divisions.motherInd] == idx);
    for j=1:length(dd)
                        
        ff = divisions(dd(j)).motherFrame;        
        divMatrix(idx,ff)=1;
        
    end
        
    dd = find([divisions.sisterInd] == idx);    
    for j=1:length(dd)
                        
        ff = divisions(dd(j)).sisterFrame;        
        divMatrix(idx,ff)=1;        
    end
    
    %detect from trace with low false neg (hopefully)
    maxp = getDivisionsFromTrace(tmp,1,0.1);
    for j=1:size(maxp,1)
        
        %check if there is no division already close by          
        sel = (maxp(j,1)-4):(maxp(j,1)+4);
        sel = sel(sel>0); sel = sel(sel <= size(ind,2));
        
        if( sum( divMatrix(idx,sel)) == 0) 
            divMatrix(idx,maxp(j,1))=1;
        end
        
    end
    
    %correct divisions that are too close in time    
    divs = find(divMatrix(idx,:));    
    for j=1:length(divs)-1
        
        div2divTime = expe.t(divs(j+1))-expe.t(divs(j));        
        if( div2divTime < minTimeBetweenDivisions)           
            
                tmp = fillTrace(tmp);
                tmp = imnorm(tmp);

                f = [1 1 -5 1 1];
                pe = conv(tmp,f,'same');
                pe(pe<0)=0;
                pe = smooth(smooth(pe,3)); %%allow for small time differences
            
                score1 = pe(divs(j));                
                score2 = pe(divs(j+1));
                
                if(score1 >= score2)                   
                    divMatrix(idx,divs(j+1)) = 0;
                else
                    divMatrix(idx,divs(j)) = 0;
                end
                
        end
    end
end


clf
imagesc(divMatrix(longTraces,:))
z = (1-linspace(0,1,100));
cmap = [z; z; z]';
colormap(cmap)

save divMatrix.mat divMatrix

%% make small images around each cell for the trace tool

mkdirIfNotExist('snapShots')

doDraw = 0;
  
for n=1:length(longTraces)
    
    idx = longTraces(n);
    disp(100*n / length(longTraces));
    
    i = round(traj{idx}(:,1));
    j = round(traj{idx}(:,2));

    clf;

    s = 50;
    w = 2*s+1;

    %find fist non zero index
    nonZeroIdx = find(ind(idx,:) > 0);
    nonZeroIdx = nonZeroIdx(1);

    timeInd = round(linspace(1,w,length(1:NToTrack)));
    
    for k=nonZeroIdx:NToTrack
               
        a = imread(['zStackedYFP_Data/' num2str(k) '.png']);
        m = imread(['zStackedThreshCorrected/' num2str(k) '.png']);
        
        [seli selj] = getNeiInd(i(k),j(k),s,N1,N2);

        if( ind(idx,k)~=0 )
         
            sub_a = a(selj,seli);
            
            if(ind(idx,k)~=0)
                px = Me{k}; %erase other objects
                px=px(ind(idx,k)).PixelIdxList; 
                m=double(m);
                m(px)=2;
                m = m==2;
            end
            
            m = m(selj,seli);                        
         
            %make a nice gif            
            b1 = bwmorph(m,'remove');
            sub_a(b1) = min(sub_a(:));
            
            %timeline       
            if( sum( divMatrix(idx, max(1,k-2):min(N,k+2) ) > 0) )            
                sub_a(end-2:end,:) = 100;
            else
                sub_a(end-2:end,:) = 50;
            end
            
            sub_a(end-2:end,timeInd(k))=150; 
            
            %save image as png for GUI
            imwrite(sub_a,['snapShots/' num2str(idx) '_' num2str(k) '.png']);

            sub_a =  ind2rgb((sub_a)/2.0,hot);
            [imind,cm] = rgb2ind(sub_a,256);
            
            if(doDraw)
                imagesc(sub_a); drawnow;
                %pause(0.1)
            end
            
%             if k == nonZeroIdx;
%                 imwrite(imind,cm,['snapShots/' num2str(idx) '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.15);
%             else
%                 imwrite(imind,cm,['snapShots/' num2str(idx) '.gif'],'gif','WriteMode','append','DelayTime',0.15);
%             end
%             
        else
%             sub_a = zeros(w,w);            
%             
%             sub_a =  ind2rgb((sub_a)/3,gray);
%             [imind,cm] = rgb2ind(sub_a,256);
%             if k == 1;
%                 imwrite(imind,cm,['snapShots/' num2str(idx) '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.15);
%             else
%                 imwrite(imind,cm,['snapShots/' num2str(idx) '.gif'],'gif','WriteMode','append','DelayTime',0.15);
%             end            
        end
    end
end


%% Traces tool, correct divs and peaks

edit code/guiTraces.m


%% Load corrected peaks and divs

load divMatrixFinal.mat
load hadPeaksChanged.mat

imagesc(peakMatrix-divMatrix)


%% div Time distrib

p = 0.5*getBinarySeq(peakMatrix,divMatrix,[-1 1]);

hist(p(:,1)-p(:,2),-30:1:0)

%%

p = 0.5*getBinarySeq(peakMatrix,divMatrix,[-1 1 -1]);

hist(p(:,3)-p(:,1),0:1:30)

