function [clean_img] = CP_Cluster(noisy_img,CPparams)

max_VST_msi = max(noisy_img(:));
min_VST_msi = min(noisy_img(:));
VST_msi = (noisy_img - min_VST_msi) / (max_VST_msi - min_VST_msi);    % scale to [0, 1]
sz = size(VST_msi);
bparams.block_sz = [8,8];
bparams.overlap_sz = [6 6];
bparams.block_num = floor((sz(1:2) - bparams.overlap_sz)./(bparams.block_sz - bparams.overlap_sz));% 计算块的个数
info = bparams;
noisy_blocks = ExtractBlocks( VST_msi, bparams );

nclusters = ceil(prod(bparams.block_num)/CPparams.class);

info.nclusters = nclusters;

fkmeans_opt.careful = 1;
fprintf('Matching similar blocks...\n');
info.predenoise = 0;
X = reshape(noisy_blocks, prod(bparams.block_sz)*sz(3), size(noisy_blocks, 4))'; 

idx = fkmeans(X, nclusters, fkmeans_opt); 
clean_blocks = zeros(size(noisy_blocks));

for k = 1:nclusters
    matched_blocks = tensor(noisy_blocks(:, :, :, idx==k));
    CP = parafac_als(matched_blocks,CPparams.atoms,struct('tol',CPparams.error,'maxiters',CPparams.k));
    clean_blocks(:,:,:,idx==k) = double(CP);
end

clean_img = JointBlocks(clean_blocks, bparams);
clean_img = clean_img * (max_VST_msi - min_VST_msi) + min_VST_msi;

end



