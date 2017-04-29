function makeGIF(X,filename)
if nargin < 2
    filename = ['GIFOutput',char(today('datetime')),'.gif'];
end

f = figure;
for n = 1:size(X,3)
      imshow(X(:,:,n),[]);
      drawnow
      frame = getframe(f);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end
end
end