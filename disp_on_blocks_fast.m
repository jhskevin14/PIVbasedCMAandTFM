function [x,y,dx,dy,pkh] = disp_on_blocks_fast(im1,im2,iblocksize,jblocksize,overlap,padding,window,method,th,N, Ni)

% DISP_ON_BLOCKS_FAST is the core of the PIV routine that finds the
% displacement field. It computes the displacements of IM2 from IM1.
% The size of the blocks (real resolution of the PIV routine) is determind by the 
% parameter BLOCKSIZE. OVERLAP determines the overlap between different
% blocks. The size of each block of IM1 and IM2 that will be crosscorrelated is 
% (BLOCKSIZE * BLOCKSIZE). IM1 and IM2 can be
% windowed using XWINDOW2D with the method indicated in WINDOW. TH is the
% thershold to determine the centroid of the crosscorrelation peak (out of order). The
% peak is computed of a window of size 2N+1. METHOD is the method used to
% compute the centre of mass. The process to determine the peak of the
% cross correlation function is iterative and uses a continous window shift
% (see Gui and Wereley, 2002)

% dx is a matrix of x-displacements, dy a matrix of y-displacements, and
% pkh is the height of the croscorrelation peak.

% by Xavier Trepat 03/2008
% modified Dhananjay Tambe 03/2009
% 
% Typical function call
%     [x,y,dx,dy,pkh] = disp_on_blocks_fast(  im1, ... % reference image
%         im2, ... % measurement image
%         Settings.resolution, ... % resolution of the PIV grid
%         Settings.resolution, ... % resolution of the PIV grid
%         Settings.overlap, ... % overlap between blocks
%         Settings.resolution, ... % padding for the PIV, usually 0
%         'hanning', ... % window
%         '2DCentroid', ...
%         0, ... % threshold of center of mass calculation
%         1, ... % Size of the window
%         4); %, ... % number of iterations
% 
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

% padding=16; % padding for the edges of the the files

sz = size(im1); % size of the image

im1 = xExpandMatrix(im1, 1, 1, padding, padding, padding, padding, 0); % pad
im2 = xExpandMatrix(im2, 1, 1, padding, padding, padding, padding, 0);

inci = round(iblocksize*(1-overlap)); % define increment
incj = round(jblocksize*(1-overlap)); % define increment

if (inci<1 || inci>iblocksize), error('Wrong Overlap in Correlation Algorithm'), end;
if (incj<1 || incj>jblocksize), error('Wrong Overlap in Correlation Algorithm'), end;

Size1 = size([1:inci:sz(1)-iblocksize+1],2); % DT
Size2 = size([1:incj:sz(2)-jblocksize+1],2); % DT
x =   zeros(Size1,Size2); % position in x
y =   zeros(Size1,Size2); % position in y
dx =  zeros(Size1,Size2); % displacements in x
dy =  zeros(Size1,Size2); % displacements in y
pkh = zeros(Size1,Size2); % height of the xcorr peak

% MAJOR LOOP
niterations=0;
for ki = 1:inci:sz(1)-iblocksize+1
    for kj =1:incj:sz(2)-jblocksize+1
        niterations=0; DX=inf; DY=inf; % initialize iterative process
        while(niterations<=Ni && abs(dx((ki+inci-1)/inci, (kj+incj-1)/incj ) -DX)>0.02 && abs(dy((ki+inci-1)/inci, (kj+incj-1)/incj) -DY)>0.02); % set iteration conditions
            niterations=niterations+1;
            DX=(dx((ki+inci-1)/inci, (kj+incj-1)/incj ));
            DY=(dy((ki+inci-1)/inci, (kj+incj-1)/incj ));
            
            im11 = im1(ki : ki+iblocksize+2*padding-1 , kj : kj+jblocksize+2*padding-1); % crop the block with padding
            im22 = im2(ki : ki+iblocksize+2*padding-1 , kj : kj+jblocksize+2*padding-1);
            
            im11=im11/mean2(im11);
            im22=im22/mean2(im22);
            
            if(DX || DY) %skip first iteration
            im11 = xsubpix_shift(im11,DX/2,DY/2); % subpixel shift of the image
            im22 = xsubpix_shift(im22,-DX/2,-DY/2);
            end
            
            im11 = im11(padding+1 : padding+iblocksize , padding+1 : padding+jblocksize); % crop the block
            im22 = im22(padding+1 : padding+iblocksize , padding+1 : padding+jblocksize);
            
            im11 = im11 - mean2(im11); % subtract the mean
            im22 = im22 - mean2(im22);                  
            
            if(window)   
                im11 = xwindow2D(im11,window); % multiply each block by a window
                im22 = xwindow2D(im22,window);                
            end
            
            c = xcorrf2(im22,im11); % / (std(im11(:))*std(im22(:))*blocksize^2); % compute the correlation function in Fourier Space and normalize
            c = real(c(1:end-1,1:end-1)); % resize for x_cntr 

            x((ki+inci-1)/inci, (kj+incj-1)/incj ) = kj+jblocksize/2;
            y((ki+inci-1)/inci, (kj+incj-1)/incj ) = ki+iblocksize/2;         
            dx((ki+inci-1)/inci, (kj+incj-1)/incj ) = DX + x_cntr(c,N,method,'X',th); % compute the displacement
            dy((ki+inci-1)/inci, (kj+incj-1)/incj ) = DY + x_cntr(c,N,method,'Y',th);   
            pkh((ki+inci-1)/inci, (kj+incj-1)/incj ) = max2(c); % store peak height
        
        end
    end
end


%% SUBFUNCTIONS

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = xExpandMatrix(A, DimRow, DimCol, padXTop,padXBottom,padYLeft,padYRight,padding)
% xExpandMatrix(A, DimX, DimY, padXTop,padXBottom,padYLeft,padYRight)
% D is A expanded DIMROW and DIMCOL times. Pads with PADDING at the periphery
% Xavier Trepat 2008

% Dimension of the final matrix before padding
FinalDimRow = DimRow*size(A,1);
FinalDimCol = DimCol*size(A,2);

% Expand the matrix
B=imresize(A,[FinalDimRow FinalDimCol]);

% pad top and bottom
sizeB=size(B);
padXTmat(1:padXTop,1:sizeB(2))=padding;
padXBmat(1:padXBottom,1:sizeB(2))=padding;
C=cat(1,padXTmat,B,padXBmat);

% pad left and right
sizeC=size(C);
padYLmat(1:sizeC(1),1:padYLeft)=padding;
padYRmat(1:sizeC(1),1:padYRight)=padding;
D=cat(2,padYLmat,C,padYRmat);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imh = xwindow2D(im,method)
% imh = xhanning2D(im,n)
% Multiplies a matrix IM by a 2D Hanning window of n*n points, where n is 
% the size of IM. Unlike the function hanning from Matlab, xhannig2D is not 
% symetric.
%
% Xavier Trepat 2007.

n=size(im,1);
m=size(im,2);
switch method
    case('hanning')   
        [xi,xj]=meshgrid(1:n,1:m);
        w = 0.25*(1-cos(2*pi*xi/n)).*(1-cos(2*pi*xj/m));
        w = w';
    case('cosn')
        dev=8;
        
        wn=meshgrid(0:m-1,0:n-1)/(m-1);
        wn = (1-cos(pi*wn).^dev);
        wm = (meshgrid(0:n-1,0:m-1)/(n-1))';
        wm = (1-cos(pi*wm).^dev);
        w = wn.*wm;
end

imh = im.*w;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im1=xsubpix_shift(im,Xj,Xi)

Xj=-Xj;
Xi=-Xi;

I=floor(Xi);
J=floor(Xj);

xi = Xi - I;
xj = Xj - J;

im1 = size(im); % Preallocate for speed
for ki=1:size(im,1)
    for kj=1:size(im,2)
        if (ki+I>1) && (ki+I+1<size(im,1)) && (kj+J>1) && (kj+J+1<size(im,2))
            im1(ki,kj) = (1-xi)*(1-xj)*im(ki+I,  kj+J)   + ...
                         xi*(1-xj)*    im(ki+I+1,kj+J)   + ...
                         xj*(1-xi)*    im(ki+I,  kj+J+1) + ...
                         xi*xj*        im(ki+I+1,kj+J+1);
        else
            im1(ki,kj)=im(ki,kj);
        end
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = xcorrf2(a,b,pad)
%  c = xcorrf2(a,b)
%   Two-dimensional cross-correlation using Fourier transforms.
%       XCORRF2(A,B) computes the crosscorrelation of matrices A and B.
%       XCORRF2(A) is the autocorrelation function.
%       This routine is functionally equivalent to xcorr2 but usually faster.
%       See also XCORR2.

%       Author(s): R. Johnson
%       $Revision: 1.0 $  $Date: 1995/11/27 $

  if nargin<3
    pad='yes';
  end
  
  
  [ma,na] = size(a);
  if nargin == 1
    %       for autocorrelation
    b = a;
  end
  [mb,nb] = size(b);
  %       make reverse conjugate of one array
  b = conj(b(mb:-1:1,nb:-1:1));
  
  if strcmp(pad,'yes')
    %       use power of 2 transform lengths
    mf = 2^nextpow2(ma+mb);
    nf = 2^nextpow2(na+nb);
%     at = fft2(b,mf,nf);
%     bt = fft2(a,mf,nf);

    ap=[a zeros(ma,nf-na);zeros(mf-ma,na) zeros(mf-ma,nf-na)];% zeros(mf-ma,nf)];
    bp=[b zeros(mb,nf-nb);zeros(mf-mb,nb) zeros(mf-mb,nf-nb)];
 %   bp=[b zeros(mb,nb); zeros(mb,nb) zeros(mb,nb)];
    at = fft2(ap);
    bt = fft2(bp);
    
%    ai = ifft2(at);
%    bi = ifft2(bt);
    mt = size(at,1);
    nt = size(at,2);
  elseif strcmp(pad,'no')
    at = fft2(b);
    bt = fft2(a);
  
  else
    disp('Wrong input to XCORRF2'); return
  end
  
  %       multiply transforms then inverse transform
  c = ifft2(at.*bt);
  %       make real output for real input
  if ~any(any(imag(a))) && ~any(any(imag(b)))
    c = real(c);
  end
  
  % normalize 
  %c = c / (std(a(:))*std(b(:))*sqrt(ma*na*mb*nb));
  c = c/ (std(ap(:))*std(bp(:))*nt*mt);
  
  %  trim to standard size
    if strcmp(pad,'yes')
    c(ma+mb:mf,:) = [];
    c(:,na+nb:nf) = [];
  elseif strcmp(pad,'no')
    c=fftshift(c(1:end-1,1:end-1));
    
%    c(ma+mb:mf,:) = [];
%    c(:,na+nb:nf) = [];
    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = x_cntr(im,sz,method, output,th)

% X_CNTR computes the maximum of a correlation function IM with subpix accuracy
% 2*sz+1 is the size of the window over which the center is computed
% The threshold is set automatically set so that the external frame is
% always below the thhreshold.
% Output can be X or Y.
% Xavier Trepat 2007

try
    
    [mi,mj]	= find(im == max(max(im)));
    mi = mi(1);
    mj = mj(1); 
    N = size(im,1);

    switch (method),
        case '2DCentroid';
            th = max([im(max(mi-sz-1,1):min(mi+sz+1,N), max(mj-sz-1,1))'...
                    im(max(mi-sz-1,1):min(mi+sz+1,N), min(mj+sz+1,N))'...
                    im(max(mi-sz-1,1), max(mj-sz-1,1):min(mj+sz+1,N))...
                    im(min(mi+sz+1,N), max(mj-sz-1,1):min(mj+sz+1,N))]); % Set the th automatically.
            
            im1 = im(mi-sz-1:mi+sz+1, mj-sz-1:mj+sz+1); 
            
            im1 = im1-th; 
            im1(im1<0)=0;
            grdi = meshgrid(mi-sz-1 : mi+sz+1)'; 
            grdj = meshgrid(mj-sz-1 : mj+sz+1);
            
            xc = sum(sum(im1.*grdj))/sum(sum(im1));
            yc = sum(sum(im1.*grdi))/sum(sum(im1));
            
        case '1DCentroid';
            
            im1=im;
             switch output
                case 'X'
                    th=max(im1(mi,mj-sz-1),im1(mi,mj+sz+1));   
                case 'Y'
                    th=max(im1(mi-sz-1,mj),im1(mi+sz+1,mj));   
            end
            im1 = im1-th;
            im1(im1<0)=0;
            grdi = meshgrid(mi-sz : mi+sz,1)'; 
            grdj = meshgrid(mj-sz : mj+sz,1);
            xc = sum(im1(mi,mj-sz:mj+sz).*grdj)/sum(im1(mi,mj-sz:mj+sz));
            yc = sum(im1(mi-sz:mi+sz,mj).*grdi)/sum(im1(mi-sz:mi+sz,mj));
            
        case '2DPoly';
            
            im1 = im(mi-sz+1:mi+sz-1, mj-sz+1:mj+sz-1);     
            im1(im1<0)=0;
            grdi = meshgrid(mi-sz+1 : mi+sz-1)'; 
            grdj = meshgrid(mj-sz+1 : mj+sz-1);
            P = polyfitweighted2(mi-sz+1:mi+sz-1, mj-sz+1:mj+sz-1, im1, 2, ones(size(im1)));
            
            i_vals = meshgrid(min2(grdi):0.01:max2(grdi))';
            j_vals = meshgrid(min2(grdj):0.01:max2(grdj));
            im1_vals = polyval2(P,min2(grdi):0.01:max2(grdi),min2(grdj):0.01:max2(grdj));
            
            xc = j_vals(im1_vals == max2(im1_vals));
            yc = i_vals(im1_vals == max2(im1_vals));
            
        case '1DPoly';
            
            im1x = im(mi, mj-sz:mj+sz); 
            im1y = im(mi-sz:mi+sz,mj); 
            grdi = meshgrid(mi-sz : mi+sz,1)'; 
            grdj = meshgrid(mj-sz : mj+sz,1);
            switch output
                case 'X'
                    [a,b] = polyfit(grdj,im1x,2);
                    j_vals = min(grdj):0.01:max(grdj);
                    im1_vals = polyval(a, j_vals);
                    xc = j_vals(im1_vals == max(im1_vals));       
                case 'Y'                    
                    [a,b] = polyfit(grdi,im1y,2);
                    i_vals = min(grdi):0.01:max(grdi);
                    im1_vals = polyval(a, i_vals);
                    yc = i_vals(im1_vals == max(im1_vals));    
            end
    end
    
    
    if(output=='X')
        R = xc-size(im,2)/2-1;    
    elseif(output=='Y')
        R = yc-size(im,2)/2-1;    
    else
        R = NaN;
    end
    
    
catch
    R=NaN;
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = min2(a)
b = min(min(a));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = max2(a)
b = max(max(a));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%