function [bin] = obrada_slike(im)
    im_hsv = rgb2hsv(im);

%     figure();
%     subplot(3,1,1);
%     imshow(im_hsv(:,:,1));
%     title('H komponenta');
%     subplot(3,1,2);
%     imshow(im_hsv(:,:,2));
%     title('S komponenta');
%     subplot(3,1,3);
%     imshow(im_hsv(:,:,3));
%     title('V komponenta');
% 
%     figure();
%     subplot(3,1,1);
%     histogram(im_hsv(:,:,1));
%     title('Histogram H komponente');
%     xlabel('intenzitet');
%     ylabel('broj piksela');
%     subplot(3,1,2);
%     histogram(im_hsv(:,:,2));
%     title('Histogram S komponente');
%     xlabel('intenzitet');
%     ylabel('broj piksela');
%     subplot(3,1,3);
%     histogram(im_hsv(:,:,3));
%     title('Histogram V komponente');
%     xlabel('intenzitet');
%     ylabel('broj piksela');
   

    threshold = 0.2;

    bin = imbinarize(im_hsv(:,:,1), threshold);


%     figure()
%     imshow(bin)
%     title('Binarizovana slika')


    se = strel('disk', 5); 
    bin = imdilate(bin, se);
    bin = imerode(bin, se);
    [N,M] = size(bin);
    bin = ones(N,M)-bin;



%     figure()
%     imshow(bin)
%     title('Binarizovana slika nakon diletacije, erozije i inverzije')

end

