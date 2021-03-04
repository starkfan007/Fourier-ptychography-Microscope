function LED_idx = set_threshold(setofentropy, threshold)
    arraysize = sqrt(size(setofentropy, 2));
    numofrect = [8, 16, 24, 32, 40, 48, 56];
    diff = zeros(1, size(numofrect, 2))
    LED_idx = [1];
    index = [];
    for i = 1:(arraysize-1)/2
        temp = setofentropy(:,(2*i-1)^2+1:((2*i-1)^2+numofrect(:,i)));
        diff(:,i) = max(temp) - min(temp); 
        if diff(:,i) < threshold
            LED_idx = [LED_idx (2*i-1)^2+1:((2*i-1)^2+numofrect(:,i))];
        else
            index = [index (2*i-1)^2+find(temp>=0.5)];
        end
    end
    LED_idx = [LED_idx index];
end