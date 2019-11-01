
mod method;
mod linear;

use clap::{Arg, App};
use image;

#[macro_use]
extern crate log;
extern crate simplelog;
use simplelog::*;



fn main() {
    TermLogger::init(LevelFilter::Info, Config::default(), TerminalMode::Mixed).unwrap();
    let matches = App::new("image mattiing")
                    .arg(Arg::with_name("input_original_image")
                        .short("i")
                        .long("original")
                        .value_name("ORIGINAL_IMAGE")
                        .takes_value(true))
                    .arg(Arg::with_name("input_scribble_image")
                        .short("s")
                        .long("scribble")
                        .value_name("SCRIBBLE_IMAGE")
                        .takes_value(true))
                    .get_matches();
    info!("current directory: {:?}", std::env::current_dir().unwrap());
    let input_img = matches.value_of("input_original_image").unwrap();
    let s_img = matches.value_of("input_scribble_image").unwrap();
    

    let img = image::open(input_img).unwrap().to_rgb();
    let h = img.height() as usize;
    let w = img.width() as usize;
    let rgb_mat = method::RgbMatrix::from_raw_vec(img.into_raw(), h, w);

    let s_img = image::open(s_img).unwrap().to_rgb();
    let s_rgb_mat = method::RgbMatrix::from_raw_vec(s_img.into_raw(), h, w);


    method::image_method::process(rgb_mat, s_rgb_mat, 0.000001);


}
