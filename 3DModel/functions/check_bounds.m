function bounds=check_bounds(gm)
bounds_outdoors = boundingBox(gm_outdoors);
totalX = bounds(2) - bounds(1); % maxX - minX
totalY = bounds(4) - bounds(3); % maxY - minY
totalZ = bounds(6) - bounds(5); % maxZ - minZ

fprintf('\n=== GEOMETRY BOUNDING BOX CHECK ===\n');
fprintf('X Dimensions: %.3f to %.3f meters (Total: %.3f m)\n', bounds(1), bounds(2), totalX);
fprintf('Y Dimensions: %.3f to %.3f meters (Total: %.3f m)\n', bounds(3), bounds(4), totalY);
fprintf('Z Dimensions: %.3f to %.3f meters (Total: %.3f m)\n', bounds(5), bounds(6), totalZ);
fprintf('===================================\n\n');