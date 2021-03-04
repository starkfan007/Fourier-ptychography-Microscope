function seqEntropy = calc_entropy(imseqlow)
    [~, ~, num] = size(imseqlow);
    max_value = max(max(max(imseqlow)));
    min_value = min(min(min(imseqlow)));
    seqEntropy = zeros(1, num);
    for i = 1:num
        temp = uint8(255*mat2gray(imseqlow(:,:,i), [min_value max_value]));
        seqEntropy(:,i) = entropy(temp);        
    end
end