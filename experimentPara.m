mainDir = '/Users/bieler/Desktop/matlab/23september_40deg/'; %! don't forget the / at the end

expe.Nframe = 140;
expe.Nstacks = 3;
expe.filePrefix = ['YFP'];
expe.doTrans = 0;
expe.filePrefixTrans = ['TRANS'];
expe.dt = 0.5;
expe.t = linspace(0,(expe.Nframe-1) * expe.dt,expe.Nframe);

N = expe.Nframe;

nMovies = 31;
