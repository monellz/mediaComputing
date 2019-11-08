mod lib;

use lib::*;

use clap::{Arg, App};
use image;

#[macro_use]
extern crate log;
extern crate simplelog;
use simplelog::*;


fn main() {
    TermLogger::init(LevelFilter::Info, Config::default(), TerminalMode::Mixed).unwrap();
    let matches = App::new("poisson image editing")
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
                    .arg(Arg::with_name("offset")
                        .short("offset")
                        .long("offset")
                        .value_name("OFFSET")
                        .default_value("(0,0)"))
                    .arg(Arg::with_name("type")
                        .long("type")
                        .short("type")
                        .value_name("CLONE_TYPE")
                        .takes_value(true)
                        .possible_values(&["naive", "mix_gradient"])
                        .default_value("naive"))
                    .get_matches();
    info!("current directory: {:?}", std::env::current_dir().unwrap());

    let bg = matches.value_of("background").unwrap();
    let fg = matches.value_of("foreground").unwrap();
    let mask = matches.value_of("mask").unwrap();

    let bg = image::open(bg).unwrap().to_rgb();
    let bg_h = bg.height() as usize;
    let bg_w = bg.width() as usize;
    let fg = image::open(fg).unwrap().to_rgb();
    let fg_h = fg.height() as usize;
    let fg_w = fg.width() as usize;
    let mask = image::open(mask).unwrap().to_rgb();

    let bg_mat = RgbMatrix::from_raw_vec(bg.into_raw(), bg_h, bg_w);
    let fg_mat = RgbMatrix::from_raw_vec(fg.into_raw(), fg_h, fg_w);
    let fg = Foreground::from_raw(fg_mat, mask.into_raw(), (0, 0));

    let modified_bg_mat = possion::process(bg_mat, fg, CloneType::MixGradient);

    modified_bg_mat.save_img("cloned.png");
}
