function [weights] = computChannelWeight(inputs)
%COMPUTCHANNELWEIGHT 此处显示有关此函数的摘要
%   此处显示详细说明
[~, ~, c] = size(inputs);
gradients = zeros(1,c);

for i = 1 : c
    channeli = inputs(:,:,i);
    channeli = channeli ./ max(channeli(:));
    gradient_x = abs(channeli(:,2:end) - channeli(:,1:end-1));
    gradient_y = abs(channeli(2:end,:) - channeli(1:end-1,:));
    gradients(i) = sum(sum(gradient_x)) + sum(sum(gradient_y));
end

weights = gradients ./ max(gradients);
weights(weights>1) = 1.0;
end

