use std::collections::BTreeSet;

use log::*;
use image::ImageBuffer;
use linear::V3;

pub type Rgb = V3;

const MASK: u8 = 110;

pub enum ParallelType {
    Naive,
    Thread
}


pub struct RgbMatrix {
    pub rgb:    Vec<Rgb>,
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
        let mut rgb = Vec::<Rgb>::with_capacity(nrow * ncol);

        (0..raw.len()).step_by(3).for_each(|i| {
            let elem = [raw[i] as f64 / 255.0, raw[i + 1] as f64 / 255.0, raw[i + 2] as f64 / 255.0];
            rgb.push(Rgb::from_raw(elem));
        });

        RgbMatrix {
            rgb,
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
                let res = [self.rgb[elem_idx][0] * 255.0, self.rgb[elem_idx][1] * 255.0, self.rgb[elem_idx][2] * 255.0];
                *image.get_pixel_mut(i as u32, j as u32) = image::Rgb([res[0] as u8, res[1] as u8, res[2] as u8]);
            }
        }

        image.save(fname).unwrap();
    }
}

//x, y, channel
impl std::ops::Index<(usize, usize, usize)> for RgbMatrix {
    type Output = f64;
    fn index(&self, idx: (usize, usize, usize)) -> &Self::Output {
        &self.rgb[idx.0 * self.ncol + idx.1][idx.2]
    }
}

//x, y, channel
impl std::ops::IndexMut<(usize, usize, usize)> for RgbMatrix {
    fn index_mut(&mut self, idx: (usize, usize, usize)) -> &mut f64 {
        &mut self.rgb[idx.0 * self.ncol + idx.1][idx.2]
    }
}

impl std::ops::Index<(usize, usize)> for RgbMatrix {
    type Output = Rgb;
    fn index(&self, idx: (usize, usize)) -> &Self::Output {
        &self.rgb[idx.0 * self.ncol + idx.1]
    }
}

impl std::ops::IndexMut<(usize, usize)> for RgbMatrix {
    fn index_mut(&mut self, idx: (usize, usize)) -> &mut Rgb {
        &mut self.rgb[idx.0 * self.ncol + idx.1]
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


//for rayon parallel
//from https://stackoverflow.com/questions/55939552/simultaneous-mutable-access-to-arbitrary-indices-of-a-large-vector-that-are-guar
pub fn uniq_refs<'a, 'd: 'a, T> (data: &'d mut [T], indices: &'a BTreeSet<usize>) -> impl Iterator<Item = &'d mut T> + 'a {
    let start = data.as_mut_ptr();
    let in_bounds_indices = indices.range(0..data.len());

    in_bounds_indices.map(move |&i| unsafe { &mut *start.add(i) })
}