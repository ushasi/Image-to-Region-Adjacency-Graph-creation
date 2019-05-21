function imagehist = spectralhistogram(image,log_image1,log_image2,gaborfilter1,gaborfilter2,gaborfilter3,gaborfilter4,minimum,x,indexlist)
%mat(indexlist)
mat = image(:,:,1);        histedges = minimum(1) + x(1)*[0:10];
hist1 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = image(:,:,2);        histedges = minimum(2) + x(2)*[0:10];
hist2 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = image(:,:,3);        histedges = minimum(3) + x(3)*[0:10];
hist3 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = log_image1(:,:,1);        histedges = minimum(4) + x(4)*[0:10];
hist4 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = log_image1(:,:,2);        histedges = minimum(5) + x(5)*[0:10];
hist5 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = log_image1(:,:,3);        histedges = minimum(6) + x(6)*[0:10];
hist6 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = log_image2(:,:,1);        histedges = minimum(7) + x(7)*[0:10];
hist7 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = log_image2(:,:,2);        histedges = minimum(8) + x(8)*[0:10];
hist8 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = log_image2(:,:,3);        histedges = minimum(9) + x(9)*[0:10];
hist9 = histcounts(mat(indexlist),histedges,'Normalization','probability');

mat = gaborfilter1(:,:,1);        histedges = minimum(10) + x(10)*[0:10];
hist10 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter1(:,:,2);        histedges = minimum(11) + x(11)*[0:10];
hist11 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter1(:,:,3);        histedges = minimum(12) + x(12)*[0:10];
hist12 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter2(:,:,1);        histedges = minimum(13) + x(13)*[0:10];
hist13 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter2(:,:,2);        histedges = minimum(14) + x(14)*[0:10];
hist14 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter2(:,:,3);        histedges = minimum(15) + x(15)*[0:10];
hist15 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter3(:,:,1);        histedges = minimum(16) + x(16)*[0:10];
hist16 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter3(:,:,2);        histedges = minimum(17) + x(17)*[0:10];
hist17 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter3(:,:,3);        histedges = minimum(18) + x(18)*[0:10];
hist18 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter4(:,:,1);        histedges = minimum(19) + x(19)*[0:10];
hist19 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter4(:,:,2);        histedges = minimum(20) + x(20)*[0:10];
hist20 = histcounts(mat(indexlist),histedges,'Normalization','probability');
mat = gaborfilter4(:,:,3);        histedges = minimum(21) + x(21)*[0:10];
hist21 = histcounts(mat(indexlist),histedges,'Normalization','probability');

imagehist = [hist1,hist2,hist3,hist4,hist5,hist6,hist7,hist8,hist9,hist10,hist11,hist12,hist13,hist14,hist15,hist16,hist17,hist18,hist19,hist20,hist21];
clear('hist1','hist2','hist3','hist4','hist5','hist6','hist7','hist8','hist9','hist10','hist11','hist12','hist13','hist14','hist15','hist16','hist17','hist18','hist19','hist20','hist21','histedges','mat');

end