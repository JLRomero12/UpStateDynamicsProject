function [CRASTER]  =  CONV_RASTER(RASTER,width,shape)
%convolve a filter with a raster for visualization
%function [CRASTER]  =  CONV_RASTER(RASTER,width,shape)
%RASTER: INPUT matrix, each row a cell or trial, column is time
%width: INPUT S.D. (width) of the Gaussian to be used for convolving
%shape = [gauss] or step; or square = [1 1 ; 1 1];
%or exp (decay width=tau; assumes width is negative)
%CRASTER: OUTPUT convolved RASTER (matrix), the edges are already removed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RASTER CONVOLUTION
%%% when dealing with spike, makes each on into a event with width
%%% determined by GaussWidth.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified DVB 5/9/12

if (nargin==2 | strcmp(shape,'gauss'))
   GaussWidth = width;        %in ms
   filter=normpdf(-4*GaussWidth:4*GaussWidth,0,GaussWidth);
   filter=filter./max(filter);
elseif (strcmp(shape,'step'))
   filter = [1:width*4]*0;
   filter(width:width*3)=1;
elseif (strcmp(shape,'square'))
    filter=ones(width,width);
elseif (strcmp(shape,'exp'))
   span = 8*width;
   width = abs(width);
   ExpDecay = width;        %in ms
   filter = zeros(1,span*2+1);
   x=0:span;
   filter(x+span+1) = exp(-x./width);
end

%This clipping of edges is set up as 2D because the filter can be 2D
CRASTER=conv2(RASTER(:,:),filter);
[rowfilter colfilter]=size(filter);
colfilter=round((colfilter-1)/2); 
rowfilter=round((rowfilter-1)/2); 
[row col]=size(CRASTER);
CRASTER=CRASTER(:,colfilter+1:(col-colfilter));
if (rowfilter)
   CRASTER=CRASTER(rowfilter+1:(row-rowfilter),:);
end


