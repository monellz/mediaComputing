use log::*;
use image::ImageBuffer;

const MASK: u8 = 0;

pub enum CloneType {
    Naive,
    MixGradient
}

pub struct RgbMatrix {
    rgb: [Vec<f64>; 3],
    nrow:   usize,
    ncol:   usize,
}

impl RgbMatrix {
    pub fn from_raw_vec(raw: Vec<u8>, nrow: usize, ncol: usize) -> RgbMatrix {
        let mut red = Vec::<f64>::with_capacity(nrow * ncol);
        let mut green = Vec::<f64>::with_capacity(nrow * ncol);
        let mut blue = Vec::<f64>::with_capacity(nrow * ncol);

        (0..raw.len()).step_by(3).for_each(|i| {
            red.push(raw[i] as f64 / 255.0);
            green.push(raw[i + 1] as f64 / 255.0);
            blue.push(raw[i + 2] as f64 / 255.0);
        });

        RgbMatrix {
            rgb: [red, green, blue],
            nrow,
            ncol
        }
    }

    pub fn save_img(&self, fname: &str) {
        info!("save img to {}", fname);

        let mut image = ImageBuffer::<image::Rgb<u8>, Vec<_>>::new(self.ncol as u32, self.nrow as u32);
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                let elem_idx = i * self.ncol + j;
                let res = [self.rgb[0][elem_idx] * 255.0, self.rgb[1][elem_idx] * 255.0, self.rgb[2][elem_idx] * 255.0];
                *image.get_pixel_mut(i as u32, j as u32) = image::Rgb([res[0] as u8, res[1] as u8, res[2] as u8]);
            }
        }

        image.save(fname).unwrap();
    }
}


pub struct Foreground {
    mat:        RgbMatrix,
    elem_idx:   Vec<usize>,
    mask:       Vec<bool>,
    offset:     (usize, usize),
}

impl Foreground {
    pub fn from_raw(mat: RgbMatrix, mask_raw: Vec<u8>, offset: (usize, usize)) -> Foreground {
        assert_eq!(mask_raw.len(), mat.nrow * mat.ncol * 3, "mask and forground should be of the same size");

        let mut elem_idx: Vec<usize> = vec![]; 
        let mut mask: Vec<bool> = Vec::with_capacity(mat.nrow * mat.ncol);
        (0..mask_raw.len()).step_by(3).for_each(|i| {
            if mask_raw[i] == MASK {
                elem_idx.push(i / 3);
                mask.push(true);
            } else {
                mask.push(false);
            }
        });

        Foreground {
            mat,
            elem_idx,
            mask,
            offset
        }
    }

    fn local_idx(&self, p: usize) -> (usize, usize) {
        (self.elem_idx[p] / self.mat.ncol, self.elem_idx[p] % self.mat.ncol)
    }

    fn global_idx(&self, p: usize) -> (usize, usize) {
        let local = self.local_idx(p);
        (local.0 + self.offset.0, local.1 + self.offset.1)
    }


    fn is_boundary(&self, fg_idx: (i32, i32)) -> bool {
        //the obj masked should not be on the boundary of foreground
        assert!(fg_idx.0 >= 0 && fg_idx.0 < self.mat.nrow as i32 && fg_idx.1 >= 0 && fg_idx.1 < self.mat.ncol as i32, "the testing fg_idx should not be on the boundary of foreground");

        if self.mask[fg_idx.0 as usize * self.mat.ncol + fg_idx.1 as usize] {
            return false;
        }

        let dx = [0, 0, 1, -1];
        let dy = [1, -1, 0, 0];

        for i in 0..4 {
            let idx = (fg_idx.0 + dx[i], fg_idx.1 + dy[i]);
            
            if idx.0 < 0 || idx.0 >= self.mat.nrow as i32 || idx.1 < 0 || idx.1 >= self.mat.ncol as i32 {
                unreachable!("the nbh of testing fg_idx should be in the foreground")
            }

            if self.mask[idx.0 as usize * self.mat.ncol + idx.1 as usize] {
                return true;
            }
        }

        false
    }
}


pub mod possion {
    use super::*;
    use linear;
    use linear::sparse;
    pub fn process(mut bg_mat: RgbMatrix, fg: Foreground, clone_type: CloneType) -> RgbMatrix {
        let channel = ["red", "green", "blue"];
        //solve for r/g/b
        for c in 0..3 {
            info!("setup for {} channel", channel[c]);
            let (csr_mat, rhs) = setup_for_single_channel(&bg_mat, &fg, &clone_type, c);
            let cond_rev = csr_mat.get_euclid_norm_cond_rev();

            info!("calculate for {} channel", channel[c]);
            let res = csr_mat.solve_pcg_parallel(&rhs, 1e-6, 4000, None, cond_rev);

            info!("update background for {} channel", channel[c]);
            for i in 0..fg.elem_idx.len() {
                let global_idx = fg.global_idx(i);
                let global_elem_idx = global_idx.0 * bg_mat.ncol + global_idx.1;

                bg_mat.rgb[c][global_elem_idx] = res[i];
            }
        }

        bg_mat
    }

    fn setup_for_single_channel(bg_mat: &RgbMatrix, fg: &Foreground, clone_type: &CloneType, channel: usize) -> (sparse::CSRMatrix::<f64>, Vec<f64>) {
        let mut vals: Vec<f64> = vec![];
        let mut rows: Vec<usize> = vec![];
        let mut cols: Vec<usize> = vec![];
        let mut rhs: Vec<f64> = vec![];

        let dx = [-1, 1, 0, 0];
        let dy = [0, 0, 1, -1];


        (0..fg.elem_idx.len()).for_each(|p| {
            let local = fg.local_idx(p);
            vals.push(4.0);
            rows.push(p);
            cols.push(p);

            let mut b: f64 = 0.0;

            //loop for nbh
            for i in 0..4 {
                let nbh_idx = (local.0 as i32 + dx[i], local.1 as i32 + dy[i]);
                assert!(nbh_idx.0 >= 0 && nbh_idx.0 < fg.mat.nrow as i32 && nbh_idx.1 >= 0 && nbh_idx.1 < fg.mat.ncol as i32, "the nbh of inner pixel should be in the foreground");
                let q = nbh_idx.0 as usize * fg.mat.ncol + nbh_idx.1 as usize;
                b += match clone_type {
                    CloneType::Naive => fg.mat.rgb[channel][fg.elem_idx[p as usize]] - fg.mat.rgb[channel][fg.elem_idx[q as usize]],
                    CloneType::MixGradient => {
                        let delta_f = fg.mat.rgb[channel][fg.elem_idx[p as usize]] - fg.mat.rgb[channel][fg.elem_idx[q as usize]];

                        let global_p = fg.global_idx(p as usize);
                        let global_q = fg.global_idx(q as usize);
                        let global_p = global_p.0 * bg_mat.ncol + global_p.1;
                        let global_q = global_q.0 * bg_mat.ncol + global_q.1;

                        let delta_b = bg_mat.rgb[channel][global_p] - bg_mat.rgb[channel][global_q];

                        if delta_f > delta_b {
                            delta_f
                        } else {
                            delta_b
                        }
                    }
                };

                if fg.mask[q as usize] {
                    vals.push(-1.0);
                    rows.push(p as usize);
                    cols.push(q as usize);
                } else {
                    if fg.is_boundary(nbh_idx) {
                        let global_idx_q = (nbh_idx.0 as usize + fg.offset.0, nbh_idx.1 as usize + fg.offset.1);
                        assert!(global_idx_q.0 < bg_mat.nrow && global_idx_q.1 < bg_mat.ncol);

                        b += bg_mat.rgb[channel][global_idx_q.0 * bg_mat.ncol + global_idx_q.1];
                    }
                }

            }

            rhs.push(b);
        });

        let csr_mat = sparse::CSRMatrix::new(vals, rows, cols);

        (csr_mat, rhs)
    }
}