function traction
%
load('inp_setup.mat');

% Call traction_finite
load('PIV_bd.mat');
t_num = size(dx_bd,3);

% Preallocate
dx = zeros(size(x_bd));
dy = dx;
tx = dx;
ty = dx;

for t=1:t_num

    dx_t = dx_bd(:,:,t);
    dy_t = dy_bd(:,:,t);
%     dx_t(isnan(dx_t)) = 0;
%     dy_t(isnan(dy_t)) = 0;
    
    % Correct for rotations
    [dx_t, dy_t] = rotation_correction(x_bd, y_bd, dx_t, dy_t);
    
    % Subtract off median displacements to account for rigid body motion
%     dx_t = dx_t - median(dx_t(~isnan(dx_t)));
%     dy_t = dy_t - median(dy_t(~isnan(dy_t)));
    
    % Set nan displacements to zero
%     dx_t(isnan(dx_t))=0;
%     dy_t(isnan(dy_t))=0;
    
%     % If any data is large, interpolate its value
%     umag = sqrt(dx_t.^2 + dy_t.^2);
%     idx = find(abs(umag) > 5*std(umag(:))); 
%     % For a normal random variable, there's a 3e-5 chance of getting displacements greater than 4 standard deviations from zero, but the data isn't normal, so actual probability is greater than 3e-5.
%     % u_k(idx) = interp2(x,y,u_k,x(idx),y(idx));
%     % v_k(idx) = interp2(x,y,v_k,x(idx),y(idx));
%     % Setting values to nan and interpolating gives a better result for 
%     % neighboring large values than the interpolations above
%     dx_t(idx)=nan; dy_t(idx)=nan;
%     dx_t = inpaint_nans(dx_t);
%     dy_t = inpaint_nans(dy_t);
    
    % Filter displacements
    % dx_t = smooth2a(dx_t,1); % The 1 gives a 3x3 mean smoothing
    % dy_t = smooth2a(dy_t,1);
    
    % Compute tractions. Inputs are um and Pa
    [tx_t, ty_t] = traction_finite(x_bd*pix_size,...
                                   y_bd*pix_size,...
                                   dx_t*pix_size,...
                                   dy_t*pix_size,...
                                   th,...
                                   nu,...
                                   E);
    
    % Add data to displacement and traction arrays
    dx(:,:,t) = dx_t;
    dy(:,:,t) = dy_t;
    tx(:,:,t) = tx_t;
    ty(:,:,t) = ty_t;
end
    x = x_bd;
    y = y_bd;   

% Save data
save('traction.mat','x', 'y', 'dx','dy','tx','ty');
end

function [U, V] = rotation_correction(X2,Y2,UX2,UY2)
% Correction for rotations using the Kabsch algorithm.
%
% X2, Y2, UX2, UY2 : m by n matrix of coordiates and displacements
% U, V : corrected (unrotated) displacements
%
% Written by JHK 11062013
% Updated by Jacob Notbohm, Harvard School of Public Health, February 2014

[n1X2, n2X2]  = size(X2);

x2 = reshape(X2, n1X2*n2X2, 1);
y2 = reshape(Y2, n1X2*n2X2, 1);
P2 = [x2';y2';0*x2'] ;

ux2 = reshape(UX2, n1X2*n2X2, 1);
uy2 = reshape(UY2, n1X2*n2X2, 1);
Q2 = [(x2+ux2)';(y2+uy2)';0*x2'];

%%%%%%%%%%%%
%Kabsch approach
%%%%%%%%%%%%
[RR, rr, RR1] = Kabsch( Q2 , P2);
r_theta = vrrotmat2vec(RR) ;  %find corrected Q2 best matching with P2
r_theta_angle = rad2deg(r_theta(4)); % r_theta_angle: rotation angle between Q2 and P2, rr: translation
%%%%%%%%%%%%*
uu = RR1(1:2,:)-P2(1:2,:) ;
U = reshape(uu(1,:)', n1X2, n2X2 );
V = reshape(uu(2,:)', n1X2, n2X2);  % corrected displacement (rotation, translation)
end
function B=inpaint_nans(A,method)
% INPAINT_NANS: in-paints over nans in an array
% usage: B=INPAINT_NANS(A)          % default method
% usage: B=INPAINT_NANS(A,method)   % specify method used
%
% Solves approximation to one of several pdes to
% interpolate and extrapolate holes in an array
%
% arguments (input):
%   A - nxm array with some NaNs to be filled in
%
%   method - (OPTIONAL) scalar numeric flag - specifies
%       which approach (or physical metaphor to use
%       for the interpolation.) All methods are capable
%       of extrapolation, some are better than others.
%       There are also speed differences, as well as
%       accuracy differences for smooth surfaces.
%
%       methods {0,1,2} use a simple plate metaphor.
%       method  3 uses a better plate equation,
%                 but may be much slower and uses
%                 more memory.
%       method  4 uses a spring metaphor.
%       method  5 is an 8 neighbor average, with no
%                 rationale behind it compared to the
%                 other methods. I do not recommend
%                 its use.
%
%       method == 0 --> (DEFAULT) see method 1, but
%         this method does not build as large of a
%         linear system in the case of only a few
%         NaNs in a large array.
%         Extrapolation behavior is linear.
%         
%       method == 1 --> simple approach, applies del^2
%         over the entire array, then drops those parts
%         of the array which do not have any contact with
%         NaNs. Uses a least squares approach, but it
%         does not modify known values.
%         In the case of small arrays, this method is
%         quite fast as it does very little extra work.
%         Extrapolation behavior is linear.
%         
%       method == 2 --> uses del^2, but solving a direct
%         linear system of equations for nan elements.
%         This method will be the fastest possible for
%         large systems since it uses the sparsest
%         possible system of equations. Not a least
%         squares approach, so it may be least robust
%         to noise on the boundaries of any holes.
%         This method will also be least able to
%         interpolate accurately for smooth surfaces.
%         Extrapolation behavior is linear.
%
%         Note: method 2 has problems in 1-d, so this
%         method is disabled for vector inputs.
%         
%       method == 3 --+ See method 0, but uses del^4 for
%         the interpolating operator. This may result
%         in more accurate interpolations, at some cost
%         in speed.
%         
%       method == 4 --+ Uses a spring metaphor. Assumes
%         springs (with a nominal length of zero)
%         connect each node with every neighbor
%         (horizontally, vertically and diagonally)
%         Since each node tries to be like its neighbors,
%         extrapolation is as a constant function where
%         this is consistent with the neighboring nodes.
%
%       method == 5 --+ See method 2, but use an average
%         of the 8 nearest neighbors to any element.
%         This method is NOT recommended for use.
%
%
% arguments (output):
%   B - nxm array with NaNs replaced
%
%
% Example:
%  [x,y] = meshgrid(0:.01:1);
%  z0 = exp(x+y);
%  znan = z0;
%  znan(20:50,40:70) = NaN;
%  znan(30:90,5:10) = NaN;
%  znan(70:75,40:90) = NaN;
%
%  z = inpaint_nans(znan);
%
%
% See also: griddata, interp1
%
% Author: John D'Errico
% e-mail address: woodchips@rochester.rr.com
% Release: 2
% Release date: 4/15/06


% I always need to know which elements are NaN,
% and what size the array is for any method
[n,m]=size(A);
A=A(:);
nm=n*m;
k=isnan(A(:));

% list the nodes which are known, and which will
% be interpolated
nan_list=find(k);
known_list=find(~k);

% how many nans overall
nan_count=length(nan_list);

% convert NaN indices to (r,c) form
% nan_list==find(k) are the unrolled (linear) indices
% (row,column) form
[nr,nc]=ind2sub([n,m],nan_list);

% both forms of index in one array:
% column 1 == unrolled index
% column 2 == row index
% column 3 == column index
nan_list=[nan_list,nr,nc];

% supply default method
if (nargin<2) || isempty(method)
  method = 0;
elseif ~ismember(method,0:5)
  error 'If supplied, method must be one of: {0,1,2,3,4,5}.'
end

% for different methods
switch method
 case 0
  % The same as method == 1, except only work on those
  % elements which are NaN, or at least touch a NaN.
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % really a 1-d case
    work_list = nan_list(:,1);
    work_list = unique([work_list;work_list - 1;work_list + 1]);
    work_list(work_list <= 1) = [];
    work_list(work_list >= nm) = [];
    nw = numel(work_list);
    
    u = (1:nw)';
    fda = sparse(repmat(u,1,3),bsxfun(@plus,work_list,-1:1), ...
      repmat([1 -2 1],nw,1),nw,nm);
  else
    % a 2-d case
    
    % horizontal and vertical neighbors only
    talks_to = [-1 0;0 -1;1 0;0 1];
    neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
    
    % list of all nodes we have identified
    all_list=[nan_list;neighbors_list];
    
    % generate sparse array with second partials on row
    % variable for each element in either list, but only
    % for those nodes which have a row index > 1 or < n
    L = find((all_list(:,2) > 1) & (all_list(:,2) < n));
    nl=length(L);
    if nl>0
      fda=sparse(repmat(all_list(L,1),1,3), ...
        repmat(all_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1],nl,1),nm,nm);
    else
      fda=spalloc(n*m,n*m,size(all_list,1)*5);
    end
    
    % 2nd partials on column index
    L = find((all_list(:,3) > 1) & (all_list(:,3) < m));
    nl=length(L);
    if nl>0
      fda=fda+sparse(repmat(all_list(L,1),1,3), ...
        repmat(all_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
        repmat([1 -2 1],nl,1),nm,nm);
    end
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 1
  % least squares approach with del^2. Build system
  % for every array element as an unknown, and then
  % eliminate those which are knowns.

  % Build sparse matrix approximating del^2 for
  % every element in A.
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % a 1-d case
    u = (1:(nm-2))';
    fda = sparse(repmat(u,1,3),bsxfun(@plus,u,0:2), ...
      repmat([1 -2 1],nm-2,1),nm-2,nm);
  else
    % a 2-d case
    
    % Compute finite difference for second partials
    % on row variable first
    [i,j]=ndgrid(2:(n-1),1:m);
    ind=i(:)+(j(:)-1)*n;
    np=(n-2)*m;
    fda=sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
      repmat([1 -2 1],np,1),n*m,n*m);
    
    % now second partials on column variable
    [i,j]=ndgrid(1:n,2:(m-1));
    ind=i(:)+(j(:)-1)*n;
    np=n*(m-2);
    fda=fda+sparse(repmat(ind,1,3),[ind-n,ind,ind+n], ...
      repmat([1 -2 1],np,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 2
  % Direct solve for del^2 BVP across holes

  % generate sparse array with second partials on row
  % variable for each nan element, only for those nodes
  % which have a row index > 1 or < n
  
  % is it 1-d or 2-d?
  if (m == 1) || (n == 1)
    % really just a 1-d case
    error('Method 2 has problems for vector input. Please use another method.')
    
  else
    % a 2-d case
    L = find((nan_list(:,2) > 1) & (nan_list(:,2) < n));
    nl=length(L);
    if nl>0
      fda=sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-1 0 1],nl,1), ...
        repmat([1 -2 1],nl,1),n*m,n*m);
    else
      fda=spalloc(n*m,n*m,size(nan_list,1)*5);
    end
    
    % 2nd partials on column index
    L = find((nan_list(:,3) > 1) & (nan_list(:,3) < m));
    nl=length(L);
    if nl>0
      fda=fda+sparse(repmat(nan_list(L,1),1,3), ...
        repmat(nan_list(L,1),1,3)+repmat([-n 0 n],nl,1), ...
        repmat([1 -2 1],nl,1),n*m,n*m);
    end
    
    % fix boundary conditions at extreme corners
    % of the array in case there were nans there
    if ismember(1,nan_list(:,1))
      fda(1,[1 2 n+1])=[-2 1 1];
    end
    if ismember(n,nan_list(:,1))
      fda(n,[n, n-1,n+n])=[-2 1 1];
    end
    if ismember(nm-n+1,nan_list(:,1))
      fda(nm-n+1,[nm-n+1,nm-n+2,nm-n])=[-2 1 1];
    end
    if ismember(nm,nan_list(:,1))
      fda(nm,[nm,nm-1,nm-n])=[-2 1 1];
    end
    
    % eliminate knowns
    rhs=-fda(:,known_list)*A(known_list);
    
    % and solve...
    B=A;
    k=nan_list(:,1);
    B(k)=fda(k,k)\rhs(k);
    
  end
  
 case 3
  % The same as method == 0, except uses del^4 as the
  % interpolating operator.
  
  % del^4 template of neighbors
  talks_to = [-2 0;-1 -1;-1 0;-1 1;0 -2;0 -1; ...
      0 1;0 2;1 -1;1 0;1 1;2 0];
  neighbors_list=identify_neighbors(n,m,nan_list,talks_to);
  
  % list of all nodes we have identified
  all_list=[nan_list;neighbors_list];
  
  % generate sparse array with del^4, but only
  % for those nodes which have a row & column index
  % >= 3 or <= n-2
  L = find( (all_list(:,2) >= 3) & ...
            (all_list(:,2) <= (n-2)) & ...
            (all_list(:,3) >= 3) & ...
            (all_list(:,3) <= (m-2)));
  nl=length(L);
  if nl>0
    % do the entire template at once
    fda=sparse(repmat(all_list(L,1),1,13), ...
        repmat(all_list(L,1),1,13) + ...
        repmat([-2*n,-n-1,-n,-n+1,-2,-1,0,1,2,n-1,n,n+1,2*n],nl,1), ...
        repmat([1 2 -8 2 1 -8 20 -8 1 2 -8 2 1],nl,1),nm,nm);
  else
    fda=spalloc(n*m,n*m,size(all_list,1)*5);
  end
  
  % on the boundaries, reduce the order around the edges
  L = find((((all_list(:,2) == 2) | ...
             (all_list(:,2) == (n-1))) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1))) | ...
           (((all_list(:,3) == 2) | ...
             (all_list(:,3) == (m-1))) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1))));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,5), ...
      repmat(all_list(L,1),1,5) + ...
        repmat([-n,-1,0,+1,n],nl,1), ...
      repmat([1 1 -4 1 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,2) == 1) | ...
             (all_list(:,2) == n)) & ...
            (all_list(:,3) >= 2) & ...
            (all_list(:,3) <= (m-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-n,0,n],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  L = find( ((all_list(:,3) == 1) | ...
             (all_list(:,3) == m)) & ...
            (all_list(:,2) >= 2) & ...
            (all_list(:,2) <= (n-1)));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(all_list(L,1),1,3), ...
      repmat(all_list(L,1),1,3) + ...
        repmat([-1,0,1],nl,1), ...
      repmat([1 -2 1],nl,1),nm,nm);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  k=find(any(fda(:,nan_list(:,1)),2));
  
  % and solve...
  B=A;
  B(nan_list(:,1))=fda(k,nan_list(:,1))\rhs(k);
  
 case 4
  % Spring analogy
  % interpolating operator.
  
  % list of all springs between a node and a horizontal
  % or vertical neighbor
  hv_list=[-1 -1 0;1 1 0;-n 0 -1;n 0 1];
  hv_springs=[];
  for i=1:4
    hvs=nan_list+repmat(hv_list(i,:),nan_count,1);
    k=(hvs(:,2)>=1) & (hvs(:,2)<=n) & (hvs(:,3)>=1) & (hvs(:,3)<=m);
    hv_springs=[hv_springs;[nan_list(k,1),hvs(k,1)]];
  end

  % delete replicate springs
  hv_springs=unique(sort(hv_springs,2),'rows');
  
  % build sparse matrix of connections, springs
  % connecting diagonal neighbors are weaker than
  % the horizontal and vertical springs
  nhv=size(hv_springs,1);
  springs=sparse(repmat((1:nhv)',1,2),hv_springs, ...
     repmat([1 -1],nhv,1),nhv,nm);
  
  % eliminate knowns
  rhs=-springs(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  B(nan_list(:,1))=springs(:,nan_list(:,1))\rhs;
  
 case 5
  % Average of 8 nearest neighbors
  
  % generate sparse array to average 8 nearest neighbors
  % for each nan element, be careful around edges
  fda=spalloc(n*m,n*m,size(nan_list,1)*9);
  
  % -1,-1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) > 1)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,-1
  L = find(nan_list(:,3) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,-1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) > 1));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,0
  L = find(nan_list(:,2) > 1);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,0
  L = find(nan_list(:,2) < n);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % -1,+1
  L = find((nan_list(:,2) > 1) & (nan_list(:,3) < m)); 
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n-1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % 0,+1
  L = find(nan_list(:,3) < m);
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end

  % +1,+1
  L = find((nan_list(:,2) < n) & (nan_list(:,3) < m));
  nl=length(L);
  if nl>0
    fda=fda+sparse(repmat(nan_list(L,1),1,2), ...
      repmat(nan_list(L,1),1,2)+repmat([n+1, 0],nl,1), ...
      repmat([1 -1],nl,1),n*m,n*m);
  end
  
  % eliminate knowns
  rhs=-fda(:,known_list)*A(known_list);
  
  % and solve...
  B=A;
  k=nan_list(:,1);
  B(k)=fda(k,k)\rhs(k);
  
end

% all done, make sure that B is the same shape as
% A was when we came in.
B=reshape(B,n,m);

end
function neighbors_list=identify_neighbors(n,m,nan_list,talks_to)
% identify_neighbors: identifies all the neighbors of
%   those nodes in nan_list, not including the nans
%   themselves
%
% arguments (input):
%  n,m - scalar - [n,m]=size(A), where A is the
%      array to be interpolated
%  nan_list - array - list of every nan element in A
%      nan_list(i,1) == linear index of i'th nan element
%      nan_list(i,2) == row index of i'th nan element
%      nan_list(i,3) == column index of i'th nan element
%  talks_to - px2 array - defines which nodes communicate
%      with each other, i.e., which nodes are neighbors.
%
%      talks_to(i,1) - defines the offset in the row
%                      dimension of a neighbor
%      talks_to(i,2) - defines the offset in the column
%                      dimension of a neighbor
%      
%      For example, talks_to = [-1 0;0 -1;1 0;0 1]
%      means that each node talks only to its immediate
%      neighbors horizontally and vertically.
% 
% arguments(output):
%  neighbors_list - array - list of all neighbors of
%      all the nodes in nan_list

if ~isempty(nan_list)
  % use the definition of a neighbor in talks_to
  nan_count=size(nan_list,1);
  talk_count=size(talks_to,1);
  
  nn=zeros(nan_count*talk_count,2);
  j=[1,nan_count];
  for i=1:talk_count
    nn(j(1):j(2),:)=nan_list(:,2:3) + ...
        repmat(talks_to(i,:),nan_count,1);
    j=j+nan_count;
  end
  
  % drop those nodes which fall outside the bounds of the
  % original array
  L = (nn(:,1)<1)|(nn(:,1)>n)|(nn(:,2)<1)|(nn(:,2)>m); 
  nn(L,:)=[];
  
  % form the same format 3 column array as nan_list
  neighbors_list=[sub2ind([n,m],nn(:,1),nn(:,2)),nn];
  
  % delete replicates in the neighbors list
  neighbors_list=unique(neighbors_list,'rows');
  
  % and delete those which are also in the list of NaNs.
  neighbors_list=setdiff(neighbors_list,nan_list,'rows');
  
else
  neighbors_list=[];
end
end
function [U, r, PP , p0,q0, lrms] = Kabsch(P, Q, m)
% Find the Least Root Mean Square distance
% between two sets of N points in D dimensions
% and the rigid transformation (i.e. translation and rotation)
% to employ in order to bring one set that close to the other,
% Using the Kabsch (1976) algorithm.
% Note that the points are paired, i.e. we know which point in one set
% should be compared to a given point in the other set.
%
% References:
% 1) Kabsch W. A solution for the best rotation to relate two sets of vectors. Acta Cryst A 1976;32:9223.
% 2) Kabsch W. A discussion of the solution for the best rotation to relate two sets of vectors. Acta Cryst A 1978;34:8278.
% 3) http://cnx.org/content/m11608/latest/
% 4) http://en.wikipedia.org/wiki/Kabsch_algorithm
%
% We slightly generalize, allowing weights given to the points.
% Those weights are determined a priori and do not depend on the distances.
%
% We work in the convention that points are column vectors;
% some use the convention where they are row vectors instead.
%
% Input  variables:
%  P : a D*N matrix where P(a,i) is the a-th coordinate of the i-th point
%      in the 1st representation
%  Q : a D*N matrix where Q(a,i) is the a-th coordinate of the i-th point
%      in the 2nd representation
%  m : (Optional) a row vector of length N giving the weights, i.e. m(i) is
%      the weight to be assigned to the deviation of the i-th point.
%      If not supplied, we take by default the unweighted (or equal weighted)
%      m(i) = 1/N.
%      The weights do not have to be normalized;
%      we divide by the sum to ensure sum_{i=1}^N m(i) = 1.
%      The weights must be non-negative with at least one positive entry.
% Output variables:
%  U : a proper orthogonal D*D matrix, representing the rotation
%  r : a D-dimensional column vector, representing the translation
%  lrms: the Least Root Mean Square
%
% Details:
% If p_i, q_i are the i-th point (as a D-dimensional column vector)
% in the two representations, i.e. p_i = P(:,i) etc., and for
% p_i' = U p_i + r      (' does not stand for transpose!)
% we have p_i' ~ q_i, that is,
% lrms = sqrt(sum_{i=1}^N m(i) (p_i' - q_i)^2)
% is the minimal rms when going over the possible U and r.
% (assuming the weights are already normalized).
%
sz1 = size(P) ;
sz2 = size(Q) ;
if (length(sz1) ~= 2 || length(sz2) ~= 2)
    error 'P and Q must be matrices' ;
end
if (any(sz1 ~= sz2))
    error 'P and Q must be of same size' ;
end
D = sz1(1) ;         % dimension of space
N = sz1(2) ;         % number of points
if (nargin >= 3)
    if (~isvector(m) || any(size(m) ~= [1 N]))
        error 'm must be a row vector of length N' ;
    end
    if (any(m < 0))
        error 'm must have non-negative entries' ;
    end
    msum = sum(m) ;
    if (msum == 0)
        error 'm must contain some positive entry' ;
    end
    m = m / msum ;     % normalize so that weights sum to 1
else                 % m not supplied - use default
    m = ones(1,N)/N ;
end

p0 = P*m' ;          % the centroid of P
q0 = Q*m' ;          % the centroid of Q
v1 = ones(1,N) ;     % row vector of N ones
P = P - p0*v1 ;      % translating P to center the origin
Q = Q - q0*v1 ;      % translating Q to center the origin

% C is a covariance matrix of the coordinates
% C = P*diag(m)*Q'
% but this is inefficient, involving an N*N matrix, while typically D << N.
% so we use another way to compute Pdm = P*diag(m),
% which is equivalent to, but more efficient than,
% Pdm = zeros(D,N) ;
% for i=1:N
% 	Pdm(:,i) = m(i)*P(:,i) ;
% end
Pdm = bsxfun(@times,m,P) ;
C = Pdm*Q' ;
%	C = P*Q' / N ;       % (for the non-weighted case)
[V,S,W] = svd(C) ;   % singular value decomposition
I = eye(D) ;
if (det(V*W') < 0)   % more numerically stable than using (det(C) < 0)
    I(D,D) = -1 ;
end
U = W*I*V' ;

r = q0 - U*p0 ;

Diff = U*P - Q ;     % P, Q already centered
PP = U*P + q0*v1 ;   % corrected P (for rotation and translation) best matching with Q (JHK)
%	lrms = sqrt(sum(sum(Diff.*Diff))/N) ; % (for the non-weighted case)
% to compute the lrms, we employ an efficient method, equivalent to:
% lrms = 0 ;
% for i=1:N
% 	lrms = lrms + m(i)*Diff(:,i)'*Diff(:,i) ;
% end
% lrms = sqrt(lrms) ;
lrms = sqrt(sum(sum(bsxfun(@times,m,Diff).*Diff))) ;
end
function [tzx, tzy] = traction_finite(x,y,uin,vin,h,sig,E)
% TRACTION_FINITE returns a 2D traction field from a 2D displacement field.
% 
% tzx and tzy are the tractions in x and y directions
% uin and vin are the displacement fields in x and y directions
% h is the z position at which the tzx and tzy are computed (usually the height of the gel)
% d is the size of one pixel of the PIV analysis
% E is the Young's modulus
% sig is the poisson ratio
%
% UNITS:    uin, vin, h and d should be in um
%           E in Pascal
% 
% NOTE:     x and y are only used to compute the contractile moment
%
% The program implements the solution provided by Del Alamo et al (PNAS, 2007)
% A number of typos in Del Alamo et al have been identified and fixed.
%
% Xavier Trepat and Jim Butler 09/2008
%

% Compute spacing of data points (based on spatial resolution)
d = y(2,1) - y(1,1); % DT

Nalpha = size(uin,2); 
Nbeta = size(uin,1); 
s1=1-sig; % useful factors that appear later
s2=1-2*sig;
s34=3-4*sig;

U = fft2(uin).'; % Fourier transform the displacement field
V = fft2(vin).';

for kalpha=1:Nalpha, 
  %  fprintf('.')
    for kbeta=1:Nbeta,        
        
        if kalpha <= Nalpha/2, alpha = (kalpha-1)*2*pi/(Nalpha*d);
        else alpha =(-Nalpha-1+kalpha)*2*pi/(Nalpha*d);
        end
        
        if kbeta <=  Nbeta/2, beta = (kbeta-1)*2*pi/(Nbeta*d);
        else beta =(-Nbeta-1+kbeta)*2*pi/(Nbeta*d);
        end
        
        k=sqrt(alpha^2+beta^2);
        if k==0, % DC
            Tzx(1,1)=E/(2*(1+sig))*U(1,1)/h;
            Tzy(1,1)=E/(2*(1+sig))*V(1,1)/h;
        else
            uh0=U(kalpha,kbeta); % these are the fourier coeficients of the displacement field
            vh0=V(kalpha,kbeta);
            
 
            Tzx(kalpha,kbeta) = -E*beta *cosh(k*h)/(2*(1+sig)*k*sinh(k*h))*(vh0*alpha-uh0*beta) ...
                + E*alpha/(2*(1-sig^2)*k)*((s34*cosh(k*h)^2)+s2^2+(k*h)^2)/(s34*sinh(k*h)*cosh(k*h)+k*h)*(alpha*uh0+beta*vh0);
            
            Tzy(kalpha,kbeta) = +E*alpha*cosh(k*h)/(2*(1+sig)*k*sinh(k*h))*(vh0*alpha-uh0*beta) ...
                + E*beta /(2*(1-sig^2)*k)*((s34*cosh(k*h)^2)+s2^2+(k*h)^2)/(s34*sinh(k*h)*cosh(k*h)+k*h)*(alpha*uh0+beta*vh0);
               
         
        end % if k==0;
    end %for kalpha=1:Nalpha,
end % for kbeta=1:Nbeta,

% invert
tzx=real(ifft2(Tzx).');
tzy=real(ifft2(Tzy).');
end