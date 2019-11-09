use crate::*;


pub fn process(mut img: CloningImage) -> CloningImage {
    let dis2 = |a: (usize, usize), b: (usize, usize)| {
        let t = a.0 as f64 - b.0 as f64;
        t * t
    };
    let half_tan = |a: (usize, usize), b: (usize, usize), c: (usize, usize)| {
        // <a, b, c
        let ac2 = dis2(a, c);
        let bc2 = dis2(b, c);
        let ab2 = dis2(a, b);

        let cos = (ab2 + bc2 - ac2) / (2.0 * ab2.sqrt() * bc2.sqrt());

        ((1.0 - cos) / (1.0 + cos)).sqrt()
    };


    (0..img.fg_mask_idx.len()).for_each(|i| {
        if i % 100 == 0 {
            info!("pixel {}/{}", i, img.fg_mask_idx.len());
        }
        let local_i = img.fg_mask_idx[i];
        //compute the mean-value coordinates
        let mut mvc = Vec::<f64>::with_capacity(img.fg_bound_idx.len());
        let mut sum: f64 = 0.0;
        img.fg_bound_idx.iter().enumerate().for_each(|(b, &cur)| {
            let prev = if b == 0 { img.fg_bound_idx[img.fg_bound_idx.len() - 1] } else { img.fg_bound_idx[b - 1] };
            let next = if b == img.fg_bound_idx.len() - 1 { img.fg_bound_idx[0] } else { img.fg_bound_idx[b + 1] };

            let w = (half_tan(prev, local_i, cur) + half_tan(cur, local_i, next)) / dis2(cur, local_i).sqrt();

            sum += w;
            mvc.push(w);
        });

        mvc.iter_mut().for_each(|v| {
            *v /= sum;
        });


        //evaluate the mean-value interpolant 
        let mut r: [f64; 3] = [0.0, 0.0, 0.0];
        r[0] = mvc.iter().enumerate().fold(0.0, |acc, (i, &v)| acc + v * img.bound_diff[0][i]);
        r[1] = mvc.iter().enumerate().fold(0.0, |acc, (i, &v)| acc + v * img.bound_diff[1][i]);
        r[2] = mvc.iter().enumerate().fold(0.0, |acc, (i, &v)| acc + v * img.bound_diff[2][i]);
        let global_i =  img.bg_mask_idx[i];
        img.bg_mat[(0, global_i.0, global_i.1)] = r[0] + img.fg_mat[(0, local_i.0, local_i.1)]; 
        img.bg_mat[(1, global_i.0, global_i.1)] = r[1] + img.fg_mat[(1, local_i.0, local_i.1)]; 
        img.bg_mat[(2, global_i.0, global_i.1)] = r[2] + img.fg_mat[(2, local_i.0, local_i.1)]; 
    });

    img
}


pub struct CloningImage {
    pub bg_mat:     RgbMatrix,
    fg_mat:         RgbMatrix,
    fg_mask_idx:    Vec<(usize, usize)>,
    bg_mask_idx:    Vec<(usize, usize)>,


    fg_bound_idx:   Vec<(usize, usize)>,
    bound_diff:     [Vec<f64>; 3] 
}

impl CloningImage {
    pub fn from_mat(bg_mat: RgbMatrix, fg_mat: RgbMatrix, fg_mask_mat: MaskMatrix, offset: (usize, usize)) -> CloningImage {
        assert!(fg_mat.nrow == fg_mask_mat.nrow && fg_mat.ncol == fg_mask_mat.ncol, "mask and forground should be of the same size");

        let mut fg_mask_idx: Vec<(usize, usize)> = vec![];
        let mut bg_mask_idx: Vec<(usize, usize)> = vec![];

        let mut fg_bound_idx: Vec<(usize, usize)> = vec![];
        let mut bound_diff_r: Vec<f64> = vec![];
        let mut bound_diff_g: Vec<f64> = vec![];
        let mut bound_diff_b: Vec<f64> = vec![];


        let dx = [0, 0, 1, -1];
        let dy = [1, -1, 0, 0];

        for i in 0..fg_mask_mat.elem.len() {
            let local = (i / fg_mask_mat.ncol, i % fg_mask_mat.ncol);
            match fg_mask_mat.elem[i] {
                Some(_) => {
                    bg_mask_idx.push((local.0 + offset.0, local.1 + offset.1));
                    fg_mask_idx.push(local);
                },
                None => {
                    //check for boundary
                    for j in 0..4 {
                        let local_nbh = (dx[j] + local.0 as i32, dy[j] + local.1 as i32);
                        if local_nbh.0 >= 0 && local_nbh.0 < fg_mat.nrow as i32 && local_nbh.1 >= 0 && local_nbh.1 < fg_mat.ncol as i32 {
                            //nbh is in foreground
                            let local_nbh = (local_nbh.0 as usize, local_nbh.1 as usize);
                            if let Some(_) = fg_mask_mat[local_nbh] {
                                //local is in boundary
                                fg_bound_idx.push(local);

                                let global = (offset.0 + local.0 as usize, offset.1 + local.1 as usize);

                                bound_diff_r.push(bg_mat[(0, global.0, global.1)] - fg_mat[(0, local.0, local.1)]);
                                bound_diff_g.push(bg_mat[(1, global.0, global.1)] - fg_mat[(1, local.0, local.1)]);
                                bound_diff_b.push(bg_mat[(2, global.0, global.1)] - fg_mat[(2, local.0, local.1)]);

                                break;
                            }
                        }
                    }
                }
            };
 
        }


        CloningImage {
            bg_mat,
            fg_mat,
            fg_mask_idx,
            bg_mask_idx,

            fg_bound_idx,
            bound_diff: [bound_diff_r, bound_diff_g, bound_diff_b]
        }
    }
}