function [im_k, xedge, yedge, domain] = find_boundary(im_k, filter_th, dilate_rad)
    [M, N] = size(im_k);

    % Get rid of some noise in the image
    im_k = medfilt2(im_k,[3, 3]);
    
    % Find edges with sobel filter
    im_1 = edge(im_k,'sobel',filter_th);
    
    % Dilate edges
    SE = strel('disk',dilate_rad,0); im_1 = imdilate(im_1,SE);
    
    % Median filter
    %im_1 = medfilt2(im_1,[2*dilate_rad, 2*dilate_rad]);
    
    % Clear pixels near borders
    im_1 = imclearborder(im_1,4);
    
    % Get boundaries
    B = bwboundaries(im_1);
    
    % Take out largest cell in B. This is the set of boundary coordinates
    [~, max_idx] = max(cellfun('size', B, 1));
    B_rowcol = B{max_idx};
    
    % Get boundary x and y coordinates
    xedge = B_rowcol(:,2);
    yedge = B_rowcol(:,1);
    
    % Smooth boundary data
    xedge=smooth(xedge);
    yedge=smooth(yedge);
    
    % Get all pixels inside the boundary (domain)
    domain = poly2mask(xedge,yedge,M,N);
    
    % Erode domain
    domain = imerode(domain,SE);
    
    % Repeat boundary finding procedure on eroded domain
    B = bwboundaries(domain);
    [~, max_idx] = max(cellfun('size', B, 1));
    B_rowcol = flipud(B{max_idx});
    xedge = B_rowcol(:,2);
    yedge = B_rowcol(:,1);
    xedge = smooth(xedge);
    yedge = smooth(yedge);
    domain = poly2mask(xedge,yedge,M,N);    
end