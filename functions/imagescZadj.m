%% imagescZadj - function to make imagesc-plot with z-scale adjusted so midpoint==0

function imagescZadj(rfMap)

yMax = max(max(rfMap));  yMin = min(min(rfMap));
yMx = max([abs(yMax) abs(yMin)]);
imagesc(rfMap,[-yMx yMx]);

return

