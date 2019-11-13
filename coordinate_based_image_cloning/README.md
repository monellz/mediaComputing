# Coordinates for Instant Image Cloning

Rust implementation of Coordinates for Instant Image Cloning

>Farbman, Zeev, et al. "Coordinates for instant image cloning." *ACM Transactions on Graphics (TOG)*. Vol. 28. No. 3. ACM, 2009.

## Toolchain

* stable-x86_64-unknown-linux-gnu
* rustc 1.38.0 (625451e37 2019-09-23)

## Usage

```bash
coordinate-based image cloning 

USAGE:
    coordinate_based_image_cloning [OPTIONS]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -b, --background <BACKGROUND_IMAGE>    
    -f, --foreground <FOREGROUND_IMAGE>    
    -m, --mask <MASK>                      foreground and mask should be of the same size
    -o, --output <OUTPUT_PATH>              [default: cloned.png]
    -x, --x_offset <X_OFFSET>               [default: 0]
    -y, --y_offset <Y_OFFSET>               [default: 0]
```

## Result

| Background          | Foreground          | Mask                    | Result                      |
| ------------------- | ------------------- | ----------------------- | --------------------------- |
| ![bg](img/1/bg.jpg) | ![fg](img/1/fg.jpg) | ![mask](img/1/mask.jpg) | ![cloned](img/1/cloned.png) |
| ![bg](img/2/bg.jpg) | ![fg](img/2/fg.jpg) | ![mask](img/2/mask.jpg) | ![cloned](img/2/cloned.png) |

