function [dist, weight, alt_weight] = splitting_distance(pts1, pts2, spots, links)

  dist = sqrt(bsxfun(@minus,pts1(:,1),pts2(:,1).').^2 + bsxfun(@minus,pts1(:,2),pts2(:,2).').^2);
    
  ninterm = size(pts1, 1);
  nstarts = size(pts2, 1);

  frames = unique(pts1(:, end)).';
  next_intense = NaN(ninterm, 1);
  for i = frames
    if (i == length(links))
      continue;
    end

    currents = (pts1(:,end) == i);

    tmp_links = links{i+1};
    tmp_vals = NaN(max([tmp_links(:,2); pts1(currents, end-1)]), 1);
    tmp_vals(tmp_links(:,2)) = spots{i+1}(tmp_links(:,1), 3);

    next_intense(currents) = tmp_vals(pts1(currents, end-1));
  end

  weight = bsxfun(@rdivide, pts1(:,4), bsxfun(@plus,next_intense,pts2(:,4).'));
  weight(weight < 1) = weight(weight < 1).^(-2);
  weight(isnan(weight)) = Inf;
  
  %weight(isfinite(weight)) = 1;

  alt_weight = pts1(:, 4) ./ next_intense;
  alt_weight(alt_weight < 1) = alt_weight(alt_weight < 1).^(-2);
  alt_weight(isnan(alt_weight)) = max(alt_weight, [], 1);
  alt_weight = diag(alt_weight);
  %alt_weight = repmat(alt_weight, 1, ninterm);
  alt_weight(alt_weight == 0) = Inf;
  
  %alt_weight(isfinite(alt_weight)) = 100;
    
  return;
end
