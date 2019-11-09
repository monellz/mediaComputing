use log::*;
use image::ImageBuffer;


const MASK: u8 = 110;


pub struct RgbMatrix {
    pub rgb:    [Vec<f64>; 3],
    pub nrow:   usize,
    pub ncol:   usize
}

pub struct MaskMatrix {
    pub elem:   Vec<Option<usize>>,
    pub nrow:   usize,
    pub ncol:   usize
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
        let mut mask_cnt = 0 as usize;
        (0..mask.len()).step_by(3).for_each(|i| {
            if mask[i] >= MASK {
                elem.push(Some(mask_cnt));
                mask_cnt += 1;
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

