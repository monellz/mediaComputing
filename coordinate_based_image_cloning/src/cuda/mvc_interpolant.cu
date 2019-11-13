#define dis2(a_x, a_y, b_x, b_y) (((a_x) - (b_x)) * ((a_x) - (b_x)) + ((a_y) - (b_y)) * ((a_y) - (b_y)))
#define normalize(x) ((x) > 1.0? 1.0: ((x) < 0.0? 0.0: (x)))

extern "C" __device__ double half_tan(
    unsigned a_x, unsigned a_y,
    unsigned b_x, unsigned b_y,
    unsigned c_x, unsigned c_y
) {
    double ac2 = dis2((double)a_x, (double)a_y, (double)c_x, (double)c_y);
    double bc2 = dis2((double)b_x, (double)b_y, (double)c_x, (double)c_y);
    double ab2 = dis2((double)a_x, (double)a_y, (double)b_x, (double)b_y);

    double cos = (ab2 + bc2 - ac2) / (2.0 * sqrt(ab2) * sqrt(bc2));
    if (cos > 1.0) { cos = 1.0; }
    else if (cos < -1.0) { cos = -1.0; }

    return sqrt((1.0 - cos) / (1.0 + cos));
}

extern "C" __global__ void mvc_interpolant(
    //input
    const unsigned int* fg_mask_idx_x,
    const unsigned int* fg_mask_idx_y,
    const unsigned int* fg_bound_idx_x,
    const unsigned int* fg_bound_idx_y,
    const double* fg_bound_r,
    const double* fg_bound_g,
    const double* fg_bound_b,

    const unsigned int len,
    const unsigned int bound_len,

    //output
    double* out_r,
    double* out_g,
    double* out_b
) {
    //data parallel
    double sum;
    unsigned prev_x, prev_y, next_x, next_y;
    double reg_r, reg_g, reg_b;
    
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < len; i += blockDim.x * gridDim.x) {
        sum = reg_r = reg_g = reg_b = 0.0;
        for (int j = 0; j < bound_len; ++j) {
            if (j == 0) { prev_x = fg_bound_idx_x[bound_len - 1]; prev_y = fg_bound_idx_y[bound_len - 1]; }
            else { prev_x = fg_bound_idx_x[j - 1]; prev_y = fg_bound_idx_y[j - 1]; }

            if (j == bound_len - 1) { next_x = fg_bound_idx_x[0]; next_y = fg_bound_idx_y[0]; }
            else { next_x = fg_bound_idx_x[j + 1]; next_y = fg_bound_idx_y[j + 1]; }

            double w = (half_tan(prev_x, prev_y, fg_mask_idx_x[i], fg_mask_idx_y[i], fg_bound_idx_x[j], fg_bound_idx_y[j])
                        + half_tan(fg_bound_idx_x[j], fg_bound_idx_y[j], fg_mask_idx_x[i], fg_mask_idx_y[i], next_x, next_y))
                        / sqrt(dis2((double)fg_bound_idx_x[j], (double)fg_bound_idx_y[j], (double)fg_mask_idx_x[i], (double)fg_mask_idx_y[i]));

            sum += w;
            reg_r += w * fg_bound_r[j];
            reg_g += w * fg_bound_g[j];
            reg_b += w * fg_bound_b[j];
        }

        out_r[i] = reg_r / sum;
        out_g[i] = reg_g / sum;
        out_b[i] = reg_b / sum;
    }
}