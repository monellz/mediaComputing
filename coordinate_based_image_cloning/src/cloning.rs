use crate::*;
use crate::lib::uniq_refs;
use std::collections::BTreeSet;
use rayon::prelude::*;


use rustacuda::prelude::*;
use rustacuda::memory::DeviceBuffer;
use std::ffi::CString;

pub fn process(mut bg_mat: RgbMatrix, fg_mat: RgbMatrix, fg_mask_mat: MaskMatrix, offset: (usize, usize), para_type: ParallelType) -> RgbMatrix {
    let normalize = |mut x: Rgb| {
        if x[0] > 1.0 { x[0] = 1.0; }
        else if x[0] < 0.0 { x[0] = 0.0; }
        if x[1] > 1.0 { x[1] = 1.0; }
        else if x[1] < 0.0 { x[1] = 0.0; }
        if x[2] > 1.0 { x[2] = 1.0; }
        else if x[2] < 0.0 { x[2] = 0.0; }
        x
    };
    let dis2 = |a: (usize, usize), b: (usize, usize)| {
        let x = a.0 as f64 - b.0 as f64;
        let y = a.1 as f64 - b.1 as f64;
        x * x + y * y
    };
    //it can be optimized from http://pages.cs.wisc.edu/~csverma/CS777/bary.html
    let half_tan = |a: (usize, usize), b: (usize, usize), c: (usize, usize)| {
        // <a, b, c
        let ac2 = dis2(a, c);
        let bc2 = dis2(b, c);
        let ab2 = dis2(a, b);

        let cos = (ab2 + bc2 - ac2) / (2.0 * ab2.sqrt() * bc2.sqrt());
        let cos = if cos > 1.0 { 1.0f64 } else if cos < -1.0 { -1.0f64 } else { cos };

        ((1.0 - cos) / (1.0 + cos)).sqrt()
    };

    let (fg_mask_idx, fg_bound_idx, fg_bound_color) = setup_for_processing(&bg_mat, &fg_mat, &fg_mask_mat, &offset);


    match para_type {
        ParallelType::Naive => {
            info!("calculate for masked pixel, total num: {}", fg_mask_idx.len());
            fg_mask_idx.iter().enumerate().for_each(|(i, &local)| {
                if i % 1000 == 0 { info!("pixel {}/{}", i, fg_mask_idx.len()); }
                //compute the mean-value coordinates
                let mut sum: f64 = 0.0;
                let mut rgb = Rgb::new();
                fg_bound_idx.iter().enumerate().for_each(|(b, &cur)| {
                    let prev = if b == 0 { fg_bound_idx[fg_bound_idx.len() - 1] } else { fg_bound_idx[b - 1] };
                    let next = if b == fg_bound_idx.len() - 1 { fg_bound_idx[0] } else { fg_bound_idx[b + 1] };

                    let w = (half_tan(prev, local, cur) + half_tan(cur, local, next)) / dis2(cur, local).sqrt();

                    sum += w;
                    rgb += fg_bound_color[b] * w;
                });

                //evaluate the mean-value interpolant 
                let global = (local.0 + offset.0, local.1 + offset.1);
                bg_mat[global] = normalize(rgb / sum + fg_mat[local]);
            });
        },
        ParallelType::Thread => {
            info!("prepare for thread parallel");
            let mut masked_pixel: Vec<_> = {
                let mut idx = BTreeSet::new();
                fg_mask_idx.iter().for_each(|&local| {
                    let e = (local.0 + offset.0) * bg_mat.ncol + local.1 + offset.1;
                    idx.insert(e);
                });
                uniq_refs(&mut bg_mat.rgb, &idx).collect()
            };

            info!("calculate for masked pixel with threads, total num: {}", fg_mask_idx.len());
            masked_pixel.par_iter_mut().enumerate().for_each(|(i, masked_pixel)| {
                let local = fg_mask_idx[i];
                //compute the mean-value coordinates
                let mut sum: f64 = 0.0;
                let mut rgb = Rgb::new();
                fg_bound_idx.iter().enumerate().for_each(|(b, &cur)| {
                    let prev = if b == 0 { fg_bound_idx[fg_bound_idx.len() - 1] } else { fg_bound_idx[b - 1] };
                    let next = if b == fg_bound_idx.len() - 1 { fg_bound_idx[0] } else { fg_bound_idx[b + 1] };

                    let w = (half_tan(prev, local, cur) + half_tan(cur, local, next)) / dis2(cur, local).sqrt();

                    sum += w;
                    rgb += fg_bound_color[b] * w;
                });

                //evaluate the mean-value interpolant 
                **masked_pixel = normalize(rgb / sum + fg_mat[local]);
            });
        },
        ParallelType::GPU => {
            info!("prepare for GPU");
            rustacuda::init(CudaFlags::empty()).unwrap();
            let device = Device::get_device(0).unwrap();
            let _context = Context::create_and_push(ContextFlags::MAP_HOST | ContextFlags::SCHED_AUTO, device).unwrap();
            let module_data = CString::new(include_str!("./cuda/mvc_interpolant.ptx")).unwrap();
            let module = Module::load_from_string(&module_data).unwrap();
            let stream = Stream::new(StreamFlags::NON_BLOCKING, None).unwrap();

            let mut dev_fg_mask_idx_x = DeviceBuffer::from_slice(fg_mask_idx.iter().map(|x| x.0 as u32 ).collect::<Vec<_>>().as_slice()).unwrap();
            let mut dev_fg_mask_idx_y = DeviceBuffer::from_slice(fg_mask_idx.iter().map(|x| x.1 as u32 ).collect::<Vec<_>>().as_slice()).unwrap();
            let mut dev_fg_bound_idx_x = DeviceBuffer::from_slice(fg_bound_idx.iter().map(|x| x.0 as u32 ).collect::<Vec<_>>().as_slice()).unwrap();
            let mut dev_fg_bound_idx_y = DeviceBuffer::from_slice(fg_bound_idx.iter().map(|x| x.1 as u32 ).collect::<Vec<_>>().as_slice()).unwrap();
            let mut dev_fg_bound_r = DeviceBuffer::from_slice(fg_bound_color.iter().map(|x| x[0] ).collect::<Vec<_>>().as_slice()).unwrap();
            let mut dev_fg_bound_g = DeviceBuffer::from_slice(fg_bound_color.iter().map(|x| x[1] ).collect::<Vec<_>>().as_slice()).unwrap();
            let mut dev_fg_bound_b = DeviceBuffer::from_slice(fg_bound_color.iter().map(|x| x[2] ).collect::<Vec<_>>().as_slice()).unwrap();

            let mut dev_out_r = unsafe { DeviceBuffer::uninitialized(fg_mask_idx.len()).unwrap() };
            let mut dev_out_g = unsafe { DeviceBuffer::uninitialized(fg_mask_idx.len()).unwrap() };
            let mut dev_out_b = unsafe { DeviceBuffer::uninitialized(fg_mask_idx.len()).unwrap() };

            info!("call kernel function");
            unsafe {
                //launch the kernel function
                launch!(module.mvc_interpolant<<<1024, 1024, 0, stream>>>(
                    //input
                    dev_fg_mask_idx_x.as_device_ptr(),
                    dev_fg_mask_idx_y.as_device_ptr(),
                    dev_fg_bound_idx_x.as_device_ptr(),
                    dev_fg_bound_idx_y.as_device_ptr(),

                    dev_fg_bound_r.as_device_ptr(),
                    dev_fg_bound_g.as_device_ptr(),
                    dev_fg_bound_b.as_device_ptr(),

                    fg_mask_idx.len() as u32,
                    fg_bound_idx.len() as u32,

                    //output
                    dev_out_r.as_device_ptr(),
                    dev_out_g.as_device_ptr(),
                    dev_out_b.as_device_ptr()
                )).unwrap();
            }
            
            //for cpu, alloc the memory preparing for copy
            let mut out_r = vec![0.0f64; fg_mask_idx.len()];
            let mut out_g = vec![0.0f64; fg_mask_idx.len()];
            let mut out_b = vec![0.0f64; fg_mask_idx.len()];

            //also, prepare for assiging the result
            let mut masked_pixel: Vec<_> = {
                let mut idx = BTreeSet::new();
                fg_mask_idx.iter().for_each(|&local| {
                    let e = (local.0 + offset.0) * bg_mat.ncol + local.1 + offset.1;
                    idx.insert(e);
                });
                uniq_refs(&mut bg_mat.rgb[..], &idx).collect()
            };

            stream.synchronize().unwrap();

            info!("kernel synchronized");

            //copy from device to host
            dev_out_r.copy_to(&mut out_r).unwrap();
            dev_out_g.copy_to(&mut out_g).unwrap();
            dev_out_b.copy_to(&mut out_b).unwrap();

            //parallel assign
            masked_pixel.par_iter_mut().enumerate().for_each(|(i, masked_pixel)| {
               **masked_pixel = normalize(Rgb::from_raw([out_r[i], out_g[i], out_b[i]]) + fg_mat[fg_mask_idx[i]]);
            });
        }
    };

    info!("all over");
    bg_mat
}


fn setup_for_processing(bg_mat: &RgbMatrix, fg_mat: &RgbMatrix, fg_mask_mat: &MaskMatrix, offset: &(usize, usize))
    -> (Vec<(usize, usize)>, Vec<(usize, usize)>, Vec<Rgb>) {

    assert!(fg_mat.nrow == fg_mask_mat.nrow && fg_mat.ncol == fg_mask_mat.ncol, "mask and forground should be of the same size");

    info!("set up for processing");

    let mut fg_mask_idx: Vec<(usize, usize)> = vec![];

    let mut fg_bound_idx: Vec<(usize, usize)> = vec![];
    let mut fg_bound_color: Vec<Rgb> = vec![];

    let dx = [0, 1, 0, -1];
    let dy = [-1, 0, 1, 0];

    let mut boundary_map: Vec<Option<Rgb>> = vec![];

    //top point
    let mut init_boundary_pt = (fg_mask_mat.nrow, fg_mask_mat.ncol);

    let mut boundary_cnt = 0 as usize;

    for i in 0..fg_mask_mat.elem.len() {
        let local = (i / fg_mask_mat.ncol, i % fg_mask_mat.ncol);
        boundary_map.push(None);
        match fg_mask_mat.elem[i] {
            Some(_) => {
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
                            if local.0 < init_boundary_pt.0 { init_boundary_pt = local; }
                            boundary_cnt += 1;

                            let global = (offset.0 + local.0, offset.1 + local.1);

                            *boundary_map.last_mut().unwrap() = Some(bg_mat[global] - fg_mat[local]);

                            break;
                        }
                    }
                }
            }
        };
    }

    info!("bondary pixel num = {}", boundary_cnt);

    let rgb = match boundary_map[init_boundary_pt.0 * fg_mask_mat.ncol + init_boundary_pt.1] {
        Some(rgb) => rgb,
        None => unreachable!()
    };
    fg_bound_idx.push(init_boundary_pt);
    fg_bound_color.push(rgb);

    let dx = [0, 1, 1, 1, 0, -1, -1, -1];
    let dy = [-1, -1, 0, 1, 1, 1, 0, -1];

    loop {
        let mut over = true;
        for i in 0..8 {
            let nbh = (init_boundary_pt.0 as i32 + dx[i], init_boundary_pt.1 as i32 + dy[i]);
            if nbh.0 >= 0 && nbh.0 < fg_mask_mat.nrow as i32 && nbh.1 >= 0 && nbh.1 < fg_mask_mat.ncol as i32 {
                let nbh = (nbh.0 as usize, nbh.1 as usize);
                if let Some(rgb) = boundary_map[nbh.0 * fg_mask_mat.ncol + nbh.1] {
                    fg_bound_idx.push(nbh);
                    fg_bound_color.push(rgb);
                    boundary_map[init_boundary_pt.0 * fg_mask_mat.ncol + init_boundary_pt.1] = None;
                    init_boundary_pt = nbh;
                    over = false;
                    break;
                }
            }
        }
        if over { break; }
    }

    assert_eq!(boundary_cnt, fg_bound_idx.len());

    (fg_mask_idx, fg_bound_idx, fg_bound_color)
}