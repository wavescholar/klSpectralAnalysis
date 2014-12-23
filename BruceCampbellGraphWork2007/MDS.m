	[X, labels] = generate_data('helix', 2000);
	figure, scatter3(X(:,1), X(:,2), X(:,3), 5, labels); title('Original dataset'), drawnow
	no_dims = round(intrinsic_dim(X, 'MLE'));
	disp(['MLE estimate of intrinsic dimensionality: ' num2str(no_dims)]);
	mappedX = compute_mapping(X, 'Laplacian', no_dims, 7);	
	figure, scatter(mappedX(:,1), mappedX(:,2), 5, labels); title('Result of dimensionality reduction'), drawnow

    
   	[X, labels] = generate_data('3d_clusters', 2000);
	figure, scatter3(X(:,1), X(:,2), X(:,3), 5, labels); title('Original dataset'), drawnow
	no_dims = round(intrinsic_dim(X, 'MLE'));
	disp(['MLE estimate of intrinsic dimensionality: ' num2str(no_dims)]);
	mappedX = compute_mapping(X, 'Laplacian', no_dims, 5);	
	figure, scatter(mappedX(:,1), mappedX(:,2), 5, labels); title('Result of dimensionality reduction'), drawnow

    [X, labels] = generate_data('twinpeaks', 2000, 0.05);
scatter3(X(:,1), X(:,2), X(:,3), 5, labels); drawnow
d1 = intrinsic_dim(X);
disp([MLE estimation:  num2str(d1)]);
d2 = intrinsic_dim(X, 'CorrDim');
disp([Correlation dim. estimation:  num2str(d2)]);
d3 = intrinsic_dim(X, 'NearNbDim');
disp(['NN dim. estimation: ' num2str(d3)]);
d4 = intrinsic_dim(X, 'EigValue');
disp(['Eigenvalue estimation: ' num2str(d4)]);


clear all
[X, labels] = generate_data('swiss', 2000, 0.05);
d = round(intrinsic_dim(X));
Y1 = compute_mapping(X, 'LLE');
Y2 = compute_mapping(X, 'LLE', d, 7);
Y3 = compute_mapping(X, 'Laplacian', d, 7, 'JDQR');
Y4 = compute_mapping(X, 'LTSA', d, 7);
Y5 = compute_mapping(X, 'CCA', d, 'Matlab');
subplot(3, 2, 1), scatter3(X(:,1), X(:,2), X(:,3), 5, labels);
subplot(3, 2, 2), scatter(Y1(:,1), Y1(:,2), 5, labels);
subplot(3, 2, 3), scatter(Y2(:,1), Y2(:,2), 5, labels);
subplot(3, 2, 4), scatter(Y3(:,1), Y3(:,2), 5, labels);
subplot(3, 2, 5), scatter(Y4(:,1), Y4(:,2), 5, labels);
subplot(3, 2, 6), scatter(Y5(:,1), Y5(:,2), 5, labels);