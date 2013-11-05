%out = combineStack(names,N1,N2,deNoise,medianSize,weightsSegmentation,weightsData,doDraw)
function out = combineStack(names,N1,N2,Nz,deNoise,medianSize,weightsSegmentation,weightsData,doDraw)

a = zeros(N1,N2,Nz);

H = fspecial('gaussian',300,50);
    

for i=1:Nz

    name = names{i};

    a(:,:,i) = imnorm( double( imread(name) ) );

    tmp = a(:,:,i);

    q = quantile(tmp(:),0.90);
    tmp(tmp>q)=q;

    %highpass filter to homogenize the backgroud
    tmp = a(:,:,i) - imfilter(tmp,H,'symmetric'); 
    
    switch deNoise
        case 'BM3D'                    
            [NA, tmp] = BM3D(1, imnorm(tmp), 6);                 
        case 'median'                    
            tmp = medfilt2(tmp, medianSize*[1 1],'symmetric');

        case 'localNorm'

            tmp = imnorm(tmp);

            s = 10;
            b = 0*tmp;

            for p=1:1:(N1)


                seli = max(1,(p-s)):min(N1,(p+s));

                for j=1:1:(N2)

                    selj = max(1,(j-s)):min(N2,(j+s));

                    sub = tmp(seli,selj);
                    sub = sub - min(sub(:));
                    sub   = sub ./ max(max(sub(:)),0.35);

                    b(seli,selj) = b(seli,selj) + sub;

                end

                if ( doDraw && mod(p,1)==0 )
                    imagesc(b);
                    drawnow;
                end

            end

            %tmp = medfilt2(b, medianSize*[1 1]);

    end

    a(:,:,i) =  (imnorm(tmp));
            

end
        

out = zeros(size(a(:,:,1)));

if Nz>1
    for i=1:Nz

        out = out + a(:,:,i)*weightsSegmentation(i);

    end
else

    out = a(:,:,1);
end


out = imnorm(out);

if( ~strcmp( deNoise, 'localNorm') )
    out(out>0.5)=0.5;
    %out(out<quantile(out(:),0.2))=0.0;
    out = imnorm(out);
end

if doDraw    
    imagesc(out.^1)
    pause(0.1)    
end


