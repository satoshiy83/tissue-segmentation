
function result = ss_adjusted_Rand_index(partition,qartition)

partition = double(partition);
qartition = double(qartition);

n = size(partition,1);

table = partition' * qartition;
array = sum(table,2);
brray = sum(table,1);

index = 0;
for i = 1:size(table,1)
    for j = 1:size(table,2)
        if table(i,j) > 1
            index = index + nchoosek(table(i,j),2);
        end
    end
end

a = 0;
for i = 1:length(array)
    if array(i) > 1
        a = a + nchoosek(array(i),2);
    end
end
b = 0;
for j = 1:length(brray)
    if brray(j) > 1
        b = b + nchoosek(brray(j),2);
    end
end
expected_index = a * b / nchoosek(n,2);

max_index = (a + b) / 2;

result = (index - expected_index) / (max_index - expected_index);
end
