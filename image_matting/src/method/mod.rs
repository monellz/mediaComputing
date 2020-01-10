pub mod image_method;

pub const WIN_LEN: usize = 3;
pub const WIN_SIZE: usize = WIN_LEN * WIN_LEN;


use crate::linear;
use linear::{M33, V3};
use crate::{matrix_able_ops_for_whole, matrix_math_ops_for_single_elem, matrix_double_idx};

use image;
use image::ImageBuffer;
use image::RgbaImage;

pub type Rgb = V3;
pub type RgbMatrix = linear::DMatrix<Rgb>;
pub type Matrix = linear::DMatrix<f64>;

pub type ItemMap = std::collections::HashMap<(usize, usize), f64>;


matrix_able_ops_for_whole!{RgbMatrix, Rgb}
matrix_able_ops_for_whole!{Matrix, f64}

matrix_math_ops_for_single_elem!{RgbMatrix, f64}
matrix_math_ops_for_single_elem!{Matrix, f64}

matrix_double_idx!{RgbMatrix, Rgb}
matrix_double_idx!{Matrix, f64}



impl RgbMatrix {
    pub fn from_raw_vec(raw: Vec<u8>, nrow: usize, ncol: usize) -> RgbMatrix {
        let elem = (0..raw.len()).step_by(3).map(|i| {
            let mut c = Rgb::new();
            c[0] = raw[i] as f64 / 255.0;
            c[1] = raw[i + 1] as f64 / 255.0;
            c[2] = raw[i + 2] as f64 / 255.0;
            c
        }).collect();

        RgbMatrix::from_vec(elem, nrow, ncol)
    }

    pub fn from_rgba_raw_vec(raw: Vec<u8>, nrow: usize, ncol: usize) -> RgbMatrix {
        let elem = (0..raw.len()).step_by(4).map(|i| {
            let mut c = Rgb::new();
            c[0] = raw[i] as f64 / 255.0;
            c[1] = raw[i + 1] as f64 / 255.0;
            c[2] = raw[i + 2] as f64 / 255.0;
            c
        }).collect();

        RgbMatrix::from_vec(elem, nrow, ncol)
    }

    pub fn from_gray_vec(raw: Vec<f64>, nrow: usize, ncol: usize) -> RgbMatrix {
        let elem = raw.iter().map(|&x| {
            let mut c = Rgb::new();
            c[0] = x;
            c[1] = x;
            c[2] = x;
            c
        }).collect();
        RgbMatrix::from_vec(elem, nrow, ncol)
    }

    fn from_gray_matrix(gray_mat: Matrix) -> RgbMatrix {
        let nrow = gray_mat.nrow();
        let ncol = gray_mat.ncol();
        let v = gray_mat.into_vec();
        let elem = v.iter().map(|&x| {
            let mut c = Rgb::new();
            c[0] = x;
            c[1] = x;
            c[2] = x;
            c
        }).collect();
        RgbMatrix::from_vec(elem, nrow, ncol)
    }

    pub fn save_img(&self, fname: &str) {
        info!("save img to {}", fname);

        let mut image = ImageBuffer::<image::Rgb<u8>, Vec<_>>::new(self.ncol() as u32, self.nrow() as u32);
        for i in 0..self.ncol() {
            for j in 0..self.nrow() {
                let res = (self[(j, i)] * 255.0).copy_to_raw();
                *image.get_pixel_mut(i as u32, j as u32) = image::Rgb([res[0] as u8, res[1] as u8, res[2] as u8]);
            }
        }

        image.save(fname).unwrap();
    }

    fn into_raw_vec() -> (Vec<u8>, usize, usize) {
        unimplemented!()
    }
}


struct WindowRef<'a> {
    mat_ref:        &'a RgbMatrix,
    start_idx:      (usize, usize),
    cov_mat:        M33,
    mean_v:         Rgb,
    inner_mat:      Matrix,
    global_idx:     Vec<usize>,
}

struct WinRefMatrix<'a> {
    elem: Vec<WindowRef<'a>>,
    nrow: usize,
    ncol: usize,
}

impl<'a> WindowRef<'a> {
    fn new(mat_ref: &'a RgbMatrix, start_idx: (usize, usize)) -> WindowRef<'a> {
        let inner_mat = Matrix::new(WIN_SIZE, WIN_SIZE, 0.0);
        
        let mut global_idx = vec![0; WIN_SIZE];
        for i in 0..WIN_LEN {
            for j in 0..WIN_LEN {
                global_idx[i * WIN_LEN + j] = (start_idx.0 + i) * mat_ref.ncol() + start_idx.1 + j;
            }
        }
        WindowRef {
            mat_ref: mat_ref,
            start_idx: start_idx,
            cov_mat: M33::new(),
            mean_v: Rgb::new(),
            inner_mat: inner_mat,
            global_idx: global_idx
        }
    }
    fn calculate_distribution(&mut self) {
        //mean color vec
        for i in 0..WIN_LEN {
            for j in 0..WIN_LEN {
                self.mean_v += self.mat_ref[(self.start_idx.0 + i, self.start_idx.1 + j)];
            }
        }
        self.mean_v /= WIN_SIZE as f64;


        //covariance matrix 3 * 3
        /*
        for i in 0..WIN_LEN {
            for j in 0..WIN_LEN {
                let v = self.mat_ref[(self.start_idx.0 + i, self.start_idx.1 + j)] - self.mean_v;
                self.cov_mat += M33::from_v3_mul(&v, &v);
            }
        }
        self.cov_mat /= WIN_SIZE as f64;
        */
        for i in 0..WIN_LEN {
            for j in 0..WIN_LEN {
                let v = self.mat_ref[(self.start_idx.0 + i, self.start_idx.1 + j)];
                self.cov_mat += M33::from_v3_mul(&v, &v);
            }
        }
        self.cov_mat /= WIN_SIZE as f64;
        self.cov_mat -= M33::from_v3_mul(&self.mean_v, &self.mean_v);
    }
    fn calculate_inner_mat(&mut self, eps: f64) {
        let mut mid = M33::eye(eps / (WIN_SIZE as f64));
        mid += self.cov_mat;

        let mid_rev = mid.reverse();

        for i in 0..WIN_SIZE {
            for j in 0..WIN_SIZE {
                let g_i = self.get_global_elem_idx(i);
                assert_eq!(g_i, self.global_idx[i]);
                let g_j = self.get_global_elem_idx(j);
                assert_eq!(g_j, self.global_idx[j]);

                let mut v = mid_rev.quadratic(&(self.mat_ref.elem[g_i] - self.mean_v), &(self.mat_ref.elem[g_j] - self.mean_v));
                v = - (1.0 + v) / (WIN_SIZE as f64);

                //if g_i == g_j {
                if i == j {
                    v += 1.0;
                }

                self.inner_mat[(i, j)] = v;
            }
        }
    }

    fn get_global_elem_idx(&self, local_idx: usize) -> usize {
        let r = local_idx / WIN_LEN;
        let c = local_idx % WIN_LEN;

        (self.start_idx.0 + r) * self.mat_ref.ncol() + self.start_idx.1 + c
    }
}



impl<'a> WinRefMatrix<'a> {
    fn from_rgb_matrix(rgb_mat: &'a RgbMatrix, win_row_count: usize, win_col_count: usize) -> WinRefMatrix<'a> {
        let mut v: Vec<WindowRef<'a>> = Vec::<WindowRef<'a>>::with_capacity(win_row_count * win_col_count);

        for i in 0..win_row_count {
            for j in 0..win_col_count {
                v.push(WindowRef::<'a>::new(rgb_mat, (i, j)));
            }
        }

        WinRefMatrix::<'a> {
            elem: v,
            nrow: win_row_count,
            ncol: win_col_count,
        }
    }

    fn calculate_all_distribution(&mut self) {
        self.elem.iter_mut().for_each(|x| x.calculate_distribution());
    }

    fn calculate_all_inner_mat(&mut self, eps: f64) {
        self.elem.iter_mut().for_each(|x| x.calculate_inner_mat(eps));
    }

    fn create_tri_inds(&self) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
        let mut rows: Vec<usize> = Vec::with_capacity(WIN_SIZE * self.ncol * self.nrow);
        let mut cols: Vec<usize> = Vec::with_capacity(WIN_SIZE * self.ncol * self.nrow);
        let mut data: Vec<f64> = Vec::with_capacity(WIN_SIZE * self.ncol * self.nrow);

        self.elem.iter().for_each(|w| {
            for i in 0..WIN_SIZE {
                for j in 0..WIN_SIZE {
 
                    if i > j {
                        rows.push(w.global_idx[i]);
                        cols.push(w.global_idx[j]);
                        data.push(w.inner_mat[(i, j)]);

                        rows.push(w.global_idx[j]);
                        cols.push(w.global_idx[i]);
                        data.push(w.inner_mat[(i, j)]);
                    } else if i == j {
                        rows.push(w.global_idx[i]);
                        cols.push(w.global_idx[j]);
                        data.push(w.inner_mat[(i, j)]);
                    }
                }
            }
        });

        (rows, cols, data)
    }
}

pub fn vec_to_rgba(vec: std::vec::Vec<f64>, h: usize, w: usize) -> RgbaImage {
    let mut image = ImageBuffer::<image::Rgba<u8>, Vec<_>>::new(w as u32, h as u32);
    for i in 0..w {
        for j in 0..h {
            //*image.get_pixel_mut(i as u32, j as u32) = image::Rgba([res[0] as u8, res[1] as u8, res[2] as u8, 0]);
            let res = vec[j * w + i] * 255.0;
            *image.get_pixel_mut(i as u32, j as u32) = image::Rgba([res as u8, res as u8, res as u8, 255]);
        }
    }
    image
} 