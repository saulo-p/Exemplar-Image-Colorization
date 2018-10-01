function y = modeTies(x)
%Mode with ties marked as -1
[y, ~, ties] = mode(x, 2);

for t = 1:length(ties)
  if(length(ties{t}) > 1)
    idxs = [];
    for l = 1:length(ties{t})
      idxs = [idxs find(ties{t}(l)==x(t, :))];
      y(t) = x(t,min(idxs));
    end
  end
end

end