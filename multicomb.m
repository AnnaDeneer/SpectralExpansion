function combs = multicomb(vectors)

n = numel(vectors); %// number of vectors
combs = cell(1,n); %// pre-define to generate comma-separated list
[combs{end:-1:1}] = ndgrid(vectors{end:-1:1}); %// the reverse order in these two
%// comma-separated lists is needed to produce the rows of the result matrix in
%// lexicographical order 
combs = cat(n+1, combs{:}); %// concat the n n-dim arrays along dimension n+1
combs = reshape(combs,[],n);