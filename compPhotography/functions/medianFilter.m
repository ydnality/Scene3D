function [ outputImage ] = medianFilter( inputImage, kernelSize )

outputImage = inputImage;
imageWidth = size(inputImage, 2);
imageHeight = size(inputImage, 1);
for i = 1:imageWidth
    for j = 1:imageHeight
        outputImage(j,i) = median(inputImage(j, max(i - kernelSize, 1): min(i+kernelSize, imageWidth)));
    end
end
end

