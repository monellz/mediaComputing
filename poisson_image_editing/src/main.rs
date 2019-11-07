mod lib;

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
                    .arg(Arg::with_name("mask"))
                        .short("mk")
                        .long("mask")
                        .value_name("MASK")
                        .takes_value(true)
                    .get_matches();
    info!("current directory: {:?}", std::env::current_dir().unwrap());
}
