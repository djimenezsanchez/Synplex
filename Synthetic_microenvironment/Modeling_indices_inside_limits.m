function [Indices] = Modeling_indices_inside_limits(Im,Contxt_limit)
    Im = zeros(size(Im));
    Im(Contxt_limit:end-Contxt_limit,Contxt_limit:end-Contxt_limit) = 1; 
    Indices =find(Im==1);
end