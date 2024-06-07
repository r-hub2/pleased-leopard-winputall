functions {

  // Function to compute the transformation of beta
  vector mod_rp_fun(vector beta, int distrib_method) {

    int nb_rp = size(beta);
    vector[nb_rp] tbeta;
    if (distrib_method == 1) {
      tbeta = exp(beta);

    } else if (distrib_method == 2) {
      tbeta = beta;

    } else if (distrib_method == 3) {
      for(i in 1:nb_rp){
        tbeta[i] = fmax(0, beta[i] );
      }
    }
    return tbeta;
  }

  int symmat_size(int n) {
    int sz;
    // This calculates it iteratively because Stan gives a warning
    // with integer division.
    sz = 0;
    for (i in 1:n) {
      sz = sz + i;
    }
    return sz;
  }

  matrix vector_to_symmat(vector x, int n) {
    matrix[n, n] m;
    int k;
    k = 1;
    for (j in 1:n) {
      for (i in 1:j) {
        m[i, j] = x[k];
        if (i != j) {
          m[j, i] = m[i, j];
        }
        k = k + 1;
      }
    }
    return m;
  }

  vector symmat_to_vector(matrix x) {
    vector[symmat_size(rows(x))] v;
    int k;
    k = 1;
    // if x is m x n symmetric, then this will return
    // only parts of an m x m matrix.
    for (j in 1:rows(x)) {
      for (i in 1:j) {
        v[k] = x[i, j];
        k = k + 1;
      }
    }
    return v;
  }

}
