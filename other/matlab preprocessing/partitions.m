function plist = partitions(total_sum,candidate_set,max_count,fixed_count)
%
if (nargin<2) || isempty(candidate_set)
  candidate_set = 1:total_sum;
end

% how many candidates are there
n = length(candidate_set);

% error checks
if any(candidate_set<=0)
  error('All members of candidate_set must be > 0')
end
% candidates must be sorted in increasng order
if any(diff(candidate_set)<0)
  error('Efficiency requires that candidate_set be sorted')
end

% check for a max_count. do we supply a default?
if (nargin<3) || isempty(max_count)
  % how high do we need look?
  max_count = floor(total_sum./candidate_set);
elseif length(max_count)==1
  % if a scalar was provided, then turn it into a vector
  max_count = repmat(max_count,1,n);
end

% check for a fixed_count
if (nargin<4) || isempty(fixed_count)
  fixed_count = [];
elseif (fixed_count<0) || (fixed_count~=round(fixed_count))
  error('fixed_count must be a positive integer if supplied')
end

% check for degenerate cases
if isempty(fixed_count)
  if total_sum == 0
    plist = zeros(1,n);
    return
  elseif (n == 0)
    plist = [];
    return
  elseif (n == 1)
    % only one element in the set. can we form
    % total_sum from it as an integer multiple?
    p = total_sum/candidate_set;
    if (p==fix(p)) && (p<=max_count)
      plist = p;
    else
      plist = [];
    end
    return
 end
else
  % there was a fixed_count supplied
  if (total_sum == 0) && (fixed_count == 0)
    plist = zeros(1,n);
    return
  elseif (n == 0) || (fixed_count <= 0)
    plist = [];
    return
  elseif (n==1)
    % there must be a non-zero fixed_count, since
    % we did not trip the last test. since there
    % is only one candidate in the set, will it work?
    if (fixed_count*candidate_set) == total_sum
      plist = fixed_count;
    else
      plist = [];
    end
    return
  end
end

% largest element and work backwards
m = max_count(end);
c = candidate_set(end);
m = min([m,floor(total_sum/c),fixed_count]);

plist = zeros(0,n);
for i = 0:m
  temp = partitions(total_sum - i*c, ...
      candidate_set(1:(end-1)), ...
      max_count(1:(end-1)),fixed_count-i);
  plist = [plist;[temp,repmat(i,size(temp,1),1)]];  %#ok
end

