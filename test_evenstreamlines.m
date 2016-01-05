% test even streamlines
close all; clear all; clc
load wind
sx = 80;
sy = 35;
stepsize = 0.05;
maxIter = 10000;
dSep = 3;
XY = estream2(x(:,:,5),y(:,:,5),u(:,:,5),v(:,:,5),sx,sy,[stepsize,maxIter,dSep],'plot');
streamline(XY);


%% compare with matlab's stream2 function
figure(2); clf
[sx2,sy2] = meshgrid(80, 20:1:50);
streamline(stream2(x(:,:,5),y(:,:,5),u(:,:,5),v(:,:,5),sx2,sy2));

%%

 

