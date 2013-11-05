mkdirIfNotExist('Measures')

nObj=zeros(1,N);

for k=1:N
    disp(100*k/N);
    
    a = imread(['zStackedYFP_Data/' num2str(k) '.png']);
    b = imread([threshFolder num2str(k) '.png']);

    labeledImage = bwlabel(b, 8);     % Label each blob so we can make measurements of it

    Measurements = regionprops(labeledImage, a, {'PixelValues','MeanIntensity','Area','Centroid','PixelIdxList','PixelList','Eccentricity'});   

    name = ['Measures/' num2str(k) '.mat'];
            
    save(name,'Measurements')
        
    
    nObj(k)=length(Measurements);
end

clf;
plot(nObj)
ylabel('number of objects')