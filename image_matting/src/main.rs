#![allow(dead_code)]
mod method;
mod linear;

use clap::{Arg, App};
use image;
use image::gif;
use image::AnimationDecoder;

use std::fs::File;

#[macro_use]
extern crate log;
extern crate simplelog;
use simplelog::*;



fn main() {
    TermLogger::init(LevelFilter::Info, Config::default(), TerminalMode::Mixed).unwrap();
    let matches = App::new("image/gif mattiing")
                    .arg(Arg::with_name("original_image/gif")
                        .short("i")
                        .long("original")
                        .value_name("ORIGINAL_IMAGE/GIF")
                        .takes_value(true))
                    .arg(Arg::with_name("scribble_image/gif")
                        .short("s")
                        .long("scribble")
                        .value_name("SCRIBBLE_IMAGE/GIF")
                        .takes_value(true))
                    .arg(Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .value_name("OUTPUT_PATH")
                        .takes_value(true))
                    .arg(Arg::with_name("type")
                        .short("t")
                        .long("type")
                        .takes_value(true)
                        .possible_values(&["image", "gif"]))
                    .get_matches();
    //info!("current directory: {:?}", std::env::current_dir().unwrap());
    let input_img = matches.value_of("original_image/gif").unwrap();
    let s_img = matches.value_of("scribble_image/gif").unwrap();
    
    match matches.value_of("type").unwrap() {
        "image" => {
            let img = image::open(input_img).unwrap().to_rgb();
            let h = img.height() as usize;
            let w = img.width() as usize;
            let rgb_mat = method::RgbMatrix::from_raw_vec(img.into_raw(), h, w);

            let s_img = image::open(s_img).unwrap().to_rgb();
            let s_rgb_mat = method::RgbMatrix::from_raw_vec(s_img.into_raw(), h, w);

            let alpha = method::image_method::process(rgb_mat, s_rgb_mat, 0.000001);
            method::RgbMatrix::from_gray_vec(alpha, h, w).save_img("alpha_mat.png");
        },
        "gif" => {
            let file = File::open(input_img).unwrap();
            let decoder = gif::Decoder::new(file).unwrap();
            let frames = decoder.into_frames();
            let frames = frames.collect_frames().unwrap();


            let file = File::open(s_img).unwrap();
            let decoder = gif::Decoder::new(file).unwrap();
            let s_frames = decoder.into_frames();
            let s_frames = s_frames.collect_frames().unwrap();

            let mut res_frames = Vec::<image::Frame>::new();

            println!("frame len = {:?} s_frames = {:?}", frames.len(), s_frames.len());

            frames.into_iter().zip(s_frames.into_iter()).for_each(|(f, s_f)| {
                let left = f.left();
                let top = f.top();
                let delay = f.delay();
                
                let buf = f.into_buffer();
                let h = buf.height() as usize;
                let w = buf.width() as usize;
                let rgb_mat = method::RgbMatrix::from_rgba_raw_vec(buf.into_raw(), h, w);
                let s_rgb_mat = method::RgbMatrix::from_rgba_raw_vec(s_f.into_buffer().into_raw(), h, w);
                let alpha = method::image_method::process(rgb_mat, s_rgb_mat, 0.000001);

                let img_buf = method::vec_to_rgba(alpha, h, w);
                
                res_frames.push(image::Frame::from_parts(img_buf, left, top, delay));
            });

            let file = File::create("output.gif").unwrap();
            let mut encoder = gif::Encoder::new(file);
            encoder.encode_frames(res_frames).unwrap();
        }
        _ => {
            unreachable!("type not support");
        }
    };

}
