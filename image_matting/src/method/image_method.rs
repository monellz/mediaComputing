use sprs::{TriMat, CsMat};

use sprs_ldl::LdlNumeric;

use nalgebra::sparse::{CsMatrix, CsCholesky};
use nalgebra::sparse;

use super::*;

use crate::linear::sparse::CSRMatrix;


pub fn process(rgb_mat: RgbMatrix, s_rgb_mat: RgbMatrix, eps: f64) {
    let lambda = 100.0;


    //constrain_map: one for constrained pixel(\alpha = 0 or 1), zero for others
    let mut constrain_map = Matrix::new(rgb_mat.nrow(), rgb_mat.ncol(), 0.0);

    let mut prior = Matrix::new(rgb_mat.nrow(), rgb_mat.ncol(), 0.0);

    let diff = s_rgb_mat.clone() - rgb_mat.clone();

    for i in 0..rgb_mat.nrow() {
        for j in 0..rgb_mat.ncol() {
            let sum = diff[(i, j)].sum();
            if sum == 0.0 {
                //same color -> not constained
                constrain_map[(i, j)] = 0.0;
                prior[(i, j)] = 0.5;
            } else {
                constrain_map[(i, j)] = 1.0;
                if sum > 0.0 {
                    prior[(i, j)] = 1.0;
                } else {
                    prior[(i, j)] = 0.0;
                }
            }
        }
    }

    //(L + \lambda D_s) \alpha = \lambda b_s
    let d_s = (constrain_map.clone() * lambda).into_vec();


    let csr_mat = calculate_laplacian_csr(&rgb_mat, eps, d_s);
    //let laplacian_mat = calculate_laplacian(&rgb_mat, eps, d_s);

    //b_s: specific alpha value, one for foreground, zero for background and other pixel
    let mut b_s = Matrix::new(rgb_mat.nrow(), rgb_mat.ncol(), 0.0);

    let mut estimate = Matrix::new(rgb_mat.nrow(), rgb_mat.ncol(), 0.0);
    for i in 0..rgb_mat.nrow() {
        for j in 0..rgb_mat.ncol() {
            let sum = diff[(i, j)].sum();
            if sum > 0.0 {
                //user-applied foreground alpha = 1
                b_s[(i, j)] = 1.0 * lambda;
                estimate[(i, j)] = 1.0;
            } else if sum == 0.0 {
                estimate[(i, j)] = 0.5;
            }
        }
    }
    let estimate = estimate.into_vec();


    //now solve L \alpha = \lambda b_s
    let mut b_s = b_s.into_vec();

    /*
    let ldlt = LdlNumeric::new(laplacian_mat.view()).unwrap();
    let mut alpha = ldlt.solve(&b_s);
    */
    //let (mut alpha, mut cur_step) = csr_mat.solve_cg(&b_s, 1e-5, 100000, None);
    //let cond_rev = csr_mat.get_jacobi_cond_rev();
    //let cond_rev = csr_mat.get_abs_norm_cond_rev();
    let cond_rev = csr_mat.get_euclid_norm_cond_rev();
    //let mut alpha = csr_mat.solve_pcg(&b_s, 1e-6, 100000, None, cond_rev);
    let mut alpha = csr_mat.solve_pcg_parallel(&b_s, 1e-6, 4000, None, cond_rev);


    //nomralize
    alpha.iter_mut().for_each(|x| {
        if *x > 1.0 { *x = 1.0; }
        else if *x < 0.0 { *x = 0.0; }
    });

    //test alpha mat
    let alpha_mat = RgbMatrix::from_gray_vec(alpha, rgb_mat.nrow(), rgb_mat.ncol());
    alpha_mat.save_img("alpha_mat.png");
}

fn calculate_laplacian_csr(rgb_mat: &RgbMatrix, eps: f64, diag: Vec<f64>) -> CSRMatrix<f64> {
    let pix_num = rgb_mat.ncol() * rgb_mat.nrow();
    let win_row_count = rgb_mat.nrow() - WIN_LEN + 1;
    let win_col_count = rgb_mat.ncol() - WIN_LEN + 1;

    let mut win_mat = WinRefMatrix::from_rgb_matrix(rgb_mat, win_row_count, win_col_count);    
    
    //iterate over the windows and calculate the mean and covariance
    win_mat.calculate_all_distribution();
    info!("distribution calculated");

    win_mat.calculate_all_inner_mat(eps);
    info!("inner_mat calculated");
    
    let (mut rows, mut cols, mut data) = win_mat.create_tri_inds();
    info!("tri_inds calculated, data.len = {}", data.len());

    //add diag
    (0..diag.len()).for_each(|i| {
        rows.push(i);
        cols.push(i);
        data.push(diag[i]);
    });


    CSRMatrix::new(data, rows, cols)
}


fn calculate_laplacian(rgb_mat: &RgbMatrix, eps: f64, diag: Vec<f64>) -> CsMat<f64> {
    let pix_num = rgb_mat.ncol() * rgb_mat.nrow();
    let win_row_count = rgb_mat.nrow() - WIN_LEN + 1;
    let win_col_count = rgb_mat.ncol() - WIN_LEN + 1;
    dbg!(win_row_count * win_col_count);

    let mut win_mat = WinRefMatrix::from_rgb_matrix(rgb_mat, win_row_count, win_col_count);    
    
    //iterate over the windows and calculate the mean and covariance
    win_mat.calculate_all_distribution();
    info!("distribution calculated");

    win_mat.calculate_all_inner_mat(eps);
    info!("inner_mat calculated");

    /*
    let mut map = ItemMap::new();
    //collect all data in windows into the hashmap
    win_mat.all_reduce(&mut map);
    dbg!("all reduce over");
    */
    


    let (rows, cols, data) = win_mat.create_tri_inds();
    dbg!(data.len());
    
    let diag_rows: Vec<usize> = (0..diag.len()).collect();
    let diag_cols: Vec<usize> = (0..diag.len()).collect();
    let diag_tri_mat = TriMat::from_triplets((pix_num, pix_num), diag_rows, diag_cols, diag.to_vec()).to_csc();

    /*
    let mut data = Vec::<f64>::with_capacity(map.capacity());
    let mut rows = Vec::<usize>::with_capacity(map.capacity());
    let mut cols = Vec::<usize>::with_capacity(map.capacity());

    //add the constrain diag
    //(L + \lambda D_s) \alpha = \lambda b_s
    diag.iter().enumerate().for_each(|(i, &x)| {
        match map.get_mut(&(i, i)) {
            Some(v) => {
                *v += x;
            },
            None => {
                map.insert((i, i), x);
            },
        };
    });

    map.iter().for_each(|(&(i, j), &value)| {
        //insert all
        //TODO: all or lower?
        rows.push(i);
        cols.push(j);
        data.push(value);
        /*
        if i >= j {
            rows.push(i);
            cols.push(j);
            data.push(value);
        }
        */
    });
    */
    dbg!("tri data ready, vec len = {}", data.len());
    
    let data_tri_mat = TriMat::from_triplets((pix_num, pix_num), rows, cols, data).to_csc();

    &data_tri_mat + &diag_tri_mat
}


fn solve_sparse_system(rows: Vec<usize>, cols: Vec<usize>, data: Vec<f64>, rhs: Vec<f64>) -> Vec<f64> {
    let cs_mat = CsMatrix::from_triplet(rows.len(), cols.len(), &rows[..], &cols[..], &data[..]);
    let cholesky = CsCholesky::new(&cs_mat);
    let l = cholesky.l().unwrap();

    let rhs = nalgebra::DMatrix::from_vec(rhs.len(), 1, rhs);

    let y = l.solve_lower_triangular(&rhs).unwrap();


    let x = l.tr_solve_lower_triangular(&y).unwrap();


    x.as_slice().to_vec()
}