use log::*;
use image::ImageBuffer;

use std::collections::HashMap;

const MASK: u8 = 110;

pub enum CloneType {
    Naive,
    MixGradient
}


pub struct RgbMatrix {
    rgb:    [Vec<f64>; 3],
    nrow:   usize,
    ncol:   usize
}

pub struct MaskMatrix {
    elem:   Vec<Option<usize>>,
    nrow:   usize,
    ncol:   usize
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
        for i in 0..self.ncol {
            for j in 0..self.nrow {
                let elem_idx = j * self.ncol + i;
                let res = [self.rgb[0][elem_idx] * 255.0, self.rgb[1][elem_idx] * 255.0, self.rgb[2][elem_idx] * 255.0];
                *image.get_pixel_mut(i as u32, j as u32) = image::Rgb([res[0] as u8, res[1] as u8, res[2] as u8]);
            }
        }

        image.save(fname).unwrap();
    }
}

impl std::ops::Index<(usize, usize, usize)> for RgbMatrix {
    type Output = f64;
    fn index(&self, idx: (usize, usize, usize)) -> &Self::Output {
        &self.rgb[idx.0][idx.1 * self.ncol + idx.2]
    }
}

impl std::ops::IndexMut<(usize, usize, usize)> for RgbMatrix {
    fn index_mut(&mut self, idx: (usize, usize, usize)) -> &mut f64 {
        &mut self.rgb[idx.0][idx.1 * self.ncol + idx.2]
    }
}

impl MaskMatrix {
    pub fn from_raw_vec(mask: Vec<u8>, nrow: usize, ncol: usize) -> MaskMatrix {
        let mut elem: Vec<Option<usize>> = vec![];
        (0..mask.len()).step_by(3).for_each(|i| {
            if mask[i] >= MASK {
                elem.push(Some(elem.len()));
            } else {
                elem.push(None);
            }
        });

        MaskMatrix {
            elem,
            nrow,
            ncol
        }
    }

    pub fn save_img(&self, fname: &str) {
        info!("save mask img to {}", fname);

        let mut image = ImageBuffer::<image::Rgb<u8>, Vec<_>>::new(self.ncol as u32, self.nrow as u32);
        for i in 0..self.ncol {
            for j in 0..self.nrow {
                let elem_idx = j * self.ncol + i;
                match self.elem[elem_idx] {
                    Some(_) => *image.get_pixel_mut(i as u32, j as u32) = image::Rgb([255, 255, 255]),
                    None => *image.get_pixel_mut(i as u32, j as u32) = image::Rgb([0, 0, 0])
                };
            }
        }
        image.save(fname).unwrap();
    }
}

impl std::ops::Index<(usize, usize)> for MaskMatrix {
    type Output = Option<usize>;
    fn index(&self, idx: (usize, usize)) -> &Self::Output {
        &self.elem[idx.0 * self.ncol + idx.1]
    }
}

impl std::ops::IndexMut<(usize, usize)> for MaskMatrix {
    fn index_mut(&mut self, idx: (usize, usize)) -> &mut Option<usize> {
        &mut self.elem[idx.0 * self.ncol + idx.1]
    }
}


pub struct CloningImage {
    pub bg_mat:     RgbMatrix,
    fg_mat:         RgbMatrix,
    fg_mask_mat:    MaskMatrix,
    fg_mask_idx:    Vec<(usize, usize)>,
    bg_mask_idx:    Vec<(usize, usize)>,
    offset:         (usize, usize),
    clone_type:     CloneType,
}

impl CloningImage {
    pub fn from_mat(bg_mat: RgbMatrix, fg_mat: RgbMatrix, fg_mask_mat: MaskMatrix, offset: (usize, usize), clone_type: CloneType) -> CloningImage {
        assert!(fg_mat.nrow == fg_mask_mat.nrow && fg_mat.ncol == fg_mask_mat.ncol, "mask and forground should be of the same size");

        let mut fg_mask_idx: Vec<(usize, usize)> = vec![];
        let mut bg_mask_idx: Vec<(usize, usize)> = vec![];

        fg_mask_mat.elem.iter().enumerate().for_each(|(i, is_mask)| {
            if let Some(_) = *is_mask {
                let local = (i / fg_mask_mat.ncol, i % fg_mask_mat.ncol);
                bg_mask_idx.push((local.0 + offset.0, local.1 + offset.1));
                fg_mask_idx.push(local);
            }
        });


        CloningImage {
            bg_mat,
            fg_mat,
            fg_mask_mat,
            fg_mask_idx,
            bg_mask_idx,
            offset,
            clone_type
        }
    }


    fn is_mask_boundary(&self, global_cur: (usize, usize)) -> bool {
        //global_idx has been conformed to be in background
        
        let local_cur = (global_cur.0 as i32 - self.offset.0 as i32, global_cur.1 as i32 - self.offset.1 as i32);
        if local_cur.0 >= 0 && local_cur.1 < self.fg_mat.nrow as i32 && local_cur.1 >= 0 && local_cur.1 < self.fg_mat.ncol as i32 {
            //cur is in foreground

            //cur itself should not be masked
            if let Some(_) = self.fg_mask_mat[(local_cur.0 as usize, local_cur.1 as usize)] {
                return false;
            }
        }

        let dx = [0, 0, 1, -1];
        let dy = [1, -1, 0, 0];

        for i in 0..4 {
            let global_nbh = (dx[i] + global_cur.0 as i32, dy[i] + global_cur.1 as i32);

            if global_nbh.0 >= 0 && global_nbh.0 < self.bg_mat.nrow as i32 && global_nbh.1 >= 0 && global_nbh.1 < self.bg_mat.ncol as i32 {
                //nbh is in background
                let local_nbh = (global_nbh.0 - self.offset.0 as i32, global_nbh.1 - self.offset.1 as i32);

                if local_nbh.0 >= 0 && local_nbh.0 < self.fg_mat.nrow as i32 && local_nbh.1 >= 0 && local_nbh.1 < self.fg_mat.ncol as i32 {
                    //nbh is in foreground
                    if let Some(_) = self.fg_mask_mat[(local_nbh.0 as usize, local_nbh.1 as usize)] {
                        return true;
                    }
                }
            }
        }

        false
    } 
}



pub mod possion {
    use super::*;
    use linear;
    use linear::sparse;

    pub fn process(mut img: CloningImage) -> CloningImage {
        let channel = ["red", "green", "blue"];
        for c in 0..3 {
            info!("setup for {} channel", channel[c]);
            let (csr_mat, rhs) = setup_for_single_channel(&img, c);
            let cond_rev = csr_mat.get_euclid_norm_cond_rev();

            info!("calculate for {} channel", channel[c]);
            let res = csr_mat.solve_pcg_parallel(&rhs, 1e-6, 4000, None, cond_rev);

            info!("update background for {} channel", channel[c]);
            for i in 0..img.bg_mask_idx.len() {
                let global = img.bg_mask_idx[i];
                img.bg_mat[(c, global.0, global.1)] = res[i];
            }
        }

        img
    }

    fn setup_for_single_channel(img: &CloningImage, channel: usize) -> (sparse::CSRMatrix::<f64>, Vec<f64>) {
        let mut vals: Vec<f64> = vec![];
        let mut rows: Vec<usize> = vec![];
        let mut cols: Vec<usize> = vec![];
        let mut rhs: Vec<f64> = vec![];

        let dx = [-1, 1, 0, 0];
        let dy = [0, 0, 1, -1];

        img.fg_mask_idx.iter().enumerate().for_each(|(p, &local_p)| {
            let global_p = img.bg_mask_idx[p];

            vals.push(4.0);
            rows.push(p);
            cols.push(p);

            let mut b: f64 = 0.0;
            
            //loop for nbh
            for i in 0..4 {
                let local_nbh = (local_p.0 as i32 + dx[i], local_p.1 as i32 + dy[i]);
                let global_nbh = (global_p.0 as i32 + dx[i], global_p.1 as i32 + dy[i]);

                if local_nbh.0 < 0 || local_nbh.0 >= img.fg_mat.nrow as i32 || local_nbh.1 < 0 || local_nbh.1 >= img.fg_mat.ncol as i32 {
                    warn!("local {:?} 's nbh {:?} is not in fg_mat nrow: {} ncol: {}", local_p, local_nbh, img.fg_mat.nrow, img.fg_mat.ncol);
                }

                assert!(global_nbh.0 >= 0 && global_nbh.0 < img.bg_mat.nrow as i32 && global_nbh.1 >= 0 && global_nbh.1 < img.bg_mat.ncol as i32);

                //TODO!!
                //there are some problems!!!
                b += img.fg_mat[(channel, local_p.0, local_p.1)] - img.fg_mat[(channel, local_nbh.0 as usize, local_nbh.1 as usize)];

                //there are also some problems!!!
                if let Some(q) = img.fg_mask_mat[(local_nbh.0 as usize, local_nbh.1 as usize)] {
                    //nbh in masked field
                    vals.push(-1.0);
                    rows.push(p);
                    cols.push(q);
                } else {
                    let global_nbh = (global_nbh.0 as usize, global_nbh.1 as usize);
                    if img.is_mask_boundary(global_nbh) {
                        //nbh is in mask boundary
                        b += img.bg_mat[(channel, global_nbh.0, global_nbh.1)];
                    }
                }
            }

            rhs.push(b);
        });

        let csr_mat = sparse::CSRMatrix::new(vals, rows, cols);

        (csr_mat, rhs)
    }

}