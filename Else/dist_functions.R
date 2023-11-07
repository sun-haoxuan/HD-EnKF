dist_dim1 = function(i, j, p){return(j - i)}

dist_dim1_circ = function(i, j, p){return(min(j - i, p - j + i))}