
imDir = '/Volumes/Naef-Lab/Rosie/23september_40deg/';

experimentPara;
movies = 6;

for i=1:length(movies)
   
    n = movies(i);
    
    mkdir([mainDir 'movie' num2str(n) '/zStackedThreshCorrected/']);
    
    system(['cp -v ' imDir 'movie' num2str(n) '/zStackedThreshCorrected/*.png ' mainDir 'movie' num2str(n) '/zStackedThreshCorrected/'])
    

        
end