function img_recall(t, im_k,xedge,yedge)
sub_plot = subplot(2,1,1);
set(sub_plot, 'position', [0.025 0.05 0.95 0.95]);
imagesc(im_k); colormap gray; axis image; axis off; hold on
plot(xedge,yedge,'w');
title(['Cell image (at ' num2str((t-1)*10) ' min)']);
end