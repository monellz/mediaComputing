use clap::{Arg, App};
use image;

#[macro_use]
extern crate log;
extern crate simplelog;
use simplelog::*;


mod lib;
use lib::*;
mod cloning;


fn main() {
    TermLogger::init(LevelFilter::Info, Config::default(), TerminalMode::Mixed).unwrap();
    let matches = App::new("coordinate-based image cloning")
                    .arg(Arg::with_name("background")
                        .short("bg")
                        .long("background")
                        .value_name("BACKGROUND_IMAGE")
                        .takes_value(true))
                    .arg(Arg::with_name("foreground")
                        .short("fg")
                        .long("foreground")
                        .value_name("FOREGROUND_IMAGE")
                        .takes_value(true))
                    .arg(Arg::with_name("mask")
                        .short("mask")
                        .long("mask")
                        .value_name("MASK")
                        .takes_value(true)
                        .help("foreground and mask should be of the same size"))
                    .arg(Arg::with_name("x_offset")
                        .short("x_offset")
                        .long("x_offset")
                        .value_name("X_OFFSET")
                        .takes_value(true)
                        .default_value("0"))
                    .arg(Arg::with_name("y_offset")
                        .short("y_offset")
                        .long("y_offset")
                        .value_name("Y_OFFSET")
                        .takes_value(true)
                        .default_value("0"))
                    .arg(Arg::with_name("output")
                        .short("output")
                        .long("output")
                        .value_name("OUTPUT_PATH")
                        .takes_value(true)
                        .default_value("cloned.png"))
                    .get_matches();

    let bg = matches.value_of("background").unwrap();
    let fg = matches.value_of("foreground").unwrap();
    let mask = matches.value_of("mask").unwrap();
    let x_offset: usize = matches.value_of("x_offset").unwrap().parse().unwrap();
    let y_offset: usize = matches.value_of("y_offset").unwrap().parse().unwrap();
    let output = matches.value_of("output").unwrap();

    let bg = image::open(bg).unwrap().to_rgb();
    let bg_h = bg.height() as usize;
    let bg_w = bg.width() as usize;
    let fg = image::open(fg).unwrap().to_rgb();
    let fg_h = fg.height() as usize;
    let fg_w = fg.width() as usize;
    let mask = image::open(mask).unwrap().to_rgb();
    let mask_h = mask.height() as usize;
    let mask_w = mask.width() as usize;

    let bg_mat = RgbMatrix::from_raw_vec(bg.into_raw(), bg_h, bg_w);
    let fg_mat = RgbMatrix::from_raw_vec(fg.into_raw(), fg_h, fg_w);
    let fg_mask_mat = MaskMatrix::from_raw_vec(mask.into_raw(), mask_h, mask_w);

    let init_img = cloning::CloningImage::from_mat(bg_mat, fg_mat, fg_mask_mat, (x_offset, y_offset));

    let modified_img = cloning::process(init_img);

    modified_img.bg_mat.save_img(output);
}
