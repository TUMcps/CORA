# get image
FROM tobiasladnertum/cora:r2022a-dl

# reading ONNX networks within docker can cause as exception
# due to some gui issue (see neuralNetwork/readONNXNetwork)
# fixing it on-the-fly requires writing permission
RUN chmod +w /home/matlab/Documents/MATLAB/SupportPackages/R2022a/toolbox/nnet/supportpackages/onnx/+nnet/+internal/+cnn/+onnx/+fcn/ModelTranslation.m
RUN chmod +w /home/matlab/Documents/MATLAB/SupportPackages/R2022a/toolbox/nnet/supportpackages/onnx/+nnet/+internal/+cnn/+onnx/CustomLayerManager.m
