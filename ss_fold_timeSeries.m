
function result = ss_fold_timeSeries(data)
eata = [];
for t = 1:size(data,4)
    eata = cat(3,eata,data(:,:,:,t));
end

result = eata;
end
