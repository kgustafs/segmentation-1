mainDir = '/Users/bieler/Desktop/matlab/23september_40deg/'; %! don't forget the / at the end


imgDir = '/Volumes/Naef-Lab/Rosie/23september_40deg/';


cd(mainDir);

i = 8;

%for i=1:nMovies
   
    mkdir([imgDir 'movie' num2str(i)]);
    mkdir([imgDir 'movie' num2str(i) '/zStackedThreshSplit/']);
    mkdir([imgDir 'movie' num2str(i) '/zStackedYFP/']);
    
        
    system(['cp -v movie' num2str(i) '/links.mat '  imgDir 'movie' num2str(i) '/' ]);
    system(['cp -v movie' num2str(i) '/nuclei.mat ' imgDir 'movie' num2str(i) '/' ]);

    system(['cp -v movie' num2str(i) '/zStackedThreshSplit/*.png ' imgDir 'movie' num2str(i) '/zStackedThreshSplit/']);
    system(['cp -v movie' num2str(i) '/zStackedYFP/*.png ' imgDir 'movie' num2str(i) '/zStackedYFP/']);
        
        
%end
 