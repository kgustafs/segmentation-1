experimentPara;


imgDir = '/Volumes/Naef-Lab/Rosie/23september_40deg/';

%%
for i=1:nMovies
   
    mkdir([mainDir 'movie' num2str(i)]);
    mkdir([mainDir 'movie' num2str(i) '/img/']);
    
    if i<10
        system(['cp -v ' imgDir '00' num2str(i) '* ' mainDir 'movie' num2str(i) '/img/']);
    elseif i<100
        system(['cp -v ' imgDir '0' num2str(i) '* ' mainDir 'movie' num2str(i) '/img/']);
    else
        system(['cp -v ' imgDir '' num2str(i) '* ' mainDir 'movie' num2str(i) '/img/']);
    end
        
        
end
% 
% mkdir([mainDir 'code/']);
% system(['cp ' srcDir '*.m ' mainDir 'code/']);
% 
% %copy tracking code
% 
% mkdir([mainDir 'code/Tracking']);
% system(['cp ' '/Users/bieler/Desktop/matlab/rosie/Tracking/*.m ' mainDir 'code/Tracking']);


%% unpack tif files

doDraw = 0 ;

for i=1:nMovies

    disp(i);
    
    for j = 1:expe.Nstacks
        
        stackSuf = num2str(j);
                
        if j==1
            stackSuf='';
        end
    
        if i<10
            fname = [mainDir 'movie' num2str(i) '/img/00' num2str(i) ' YFP' stackSuf '.tif'];
        elseif i<100
            fname = [mainDir 'movie' num2str(i) '/img/0' num2str(i) ' YFP' stackSuf '.tif'];
        else
            fname = [mainDir 'movie' num2str(i) '/img/' num2str(i) ' YFP' stackSuf '.tif'];
        end

        info = imfinfo(fname);
        num_images = numel(info);
        for k = 1:num_images

            zNumber = j;
            fNumber = k;

            A = imread(fname, k);

            if(doDraw)
                imagesc(A)
                drawnow
            end

            imgName = ['YFP_' num2str(fNumber) '_' num2str(zNumber) '.TIF'];

            name = [mainDir 'movie' num2str(i) '/img/' imgName];
            %name
            imwrite(A,name)

        end
        
        system(['rm "' fname '"']);
    
    end

end
   

%%

if(expe.doTrans)

    for i=1:nMovies

        disp(i);

        if i<10
            fname = [mainDir 'movie' num2str(i) '/img/00' num2str(i) ' TRANS.tif'];
        elseif i<100
            fname = [mainDir 'movie' num2str(i) '/img/0' num2str(i) ' TRANS.tif'];
        else
            fname = [mainDir 'movie' num2str(i) '/img/' num2str(i) ' TRANS.tif'];
        end


        info = imfinfo(fname);
        num_images = numel(info);
        for k = 1:num_images
            A = imread(fname, k);

            %zNumber = mod(k-1,expe.Nstacks);
            fNumber = k;%floor((k-1)/expe.Nstacks);

            %imagesc(A)
            %drawnow

            %imgName = ['TRANS_' num2str(fNumber) '_' num2str(zNumber) '.TIF'];

            imgName = ['TRANS_' num2str(fNumber) '.TIF'];


            name = [mainDir 'movie' num2str(i) '/img/' imgName];
            %name
            imwrite(A,name)

        end

    end

end
    
