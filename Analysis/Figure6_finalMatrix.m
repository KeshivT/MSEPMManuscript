simpleMatrix = zeros(95, 93);
for a = 0:94
    for y = 1:93
        file = ['dna-methylation-changes (' num2str(y-1) ')'];
        data = readtable(file);
        m = 1;
        while (round(data{m, 1}) ~= a) & (m < height(data))
            m = m + 1;
        end
        if m ~= height(data)
            n = m;
            while round(data{n, 1}) == a
                n = n + 1;
            end
            data00 = nanmean(data(m:n, 2:3), 2);
            simpleMatrix(a + 1, y) = nanmean(table2array(data00(:, 1)));
        end
    end
end
save('simpleMatrix.mat', 'simpleMatrix');
simpMatrix = matfile('simpleMatrix.mat');
simpleMatrix = simpMatrix.simpleMatrix;
simpleMatrix(simpleMatrix == 0) = NaN;
simpleMatrix = knnimpute(simpleMatrix);
save('simpleMatrix.mat', 'simpleMatrix');