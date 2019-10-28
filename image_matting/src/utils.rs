pub const WIN_LEN: usize = 3;


pub type RgbMatrix  =   DMatrix<Rgb>;
pub type ItemMap = std::collections::HashMap<(usize, usize), f32>;


pub type Matrix     =   DMatrix<f32>;
//pub type Rgb = FixVec3;
pub type Rgb = FixRowVec3;
//pub type WinMatrix  =   DMatrix<Window>;
//    type Window     =   FixMatrix;





#[derive(Debug, Clone)]
pub struct DMatrix<T> {
    elem:   Vec<T>,
    nrow:   usize,
    ncol:   usize,
}

//this data structure is for no-copy
pub struct WindowRef<'a> {
    mat_ref:        &'a RgbMatrix,
    start_idx:      (usize, usize),
    cov_mat:        FixMatrix33,
    mean_v:         Rgb,
    inner_mat:  Matrix,
}

impl<'a> WindowRef<'a> {
    pub fn new(mat_ref: &'a RgbMatrix, start_idx: (usize, usize)) -> WindowRef<'a> {
        let lap_mat = DMatrix::<f32>::new(9, 9);
        WindowRef {
            mat_ref: mat_ref,
            start_idx: start_idx,
            cov_mat: FixMatrix33::new(),
            mean_v: Rgb::new(),
            inner_mat: lap_mat,
        }
    }

    pub fn cal_distribution(&mut self) {
        //mean color vector
        self.mean_v = Rgb::new();
        for i in 0..WIN_LEN {
            for j in 0..WIN_LEN {
                //mean_rgb = mean_rgb + self.mat_ref[(self.start_idx.0 + i, self.start_idx.1 + j)];
                //mean_rgb += self.mat_ref[(self.start_idx.0 + i, self.start_idx.1 + j)];
                self.mean_v += self.mat_ref[(self.start_idx.0 + i, self.start_idx.1 + j)];
            }
        }
        //self.mean_v = (mean_rgb / ((WIN_LEN * WIN_LEN )as f32)).into_raw();
        self.mean_v /= (WIN_LEN * WIN_LEN) as f32;
        //println!("meav_v: {:?}", self.mean_v);

        //covariance matrix 3*3
        self.cov_mat.clear();
        for i in 0..WIN_LEN {
            for j in 0..WIN_LEN {
                let v = self.mat_ref[(self.start_idx.0 + i, self.start_idx.1 + j)] - self.mean_v;
                self.cov_mat += v.transpose() * v;
            }
        }

        self.cov_mat /= (WIN_LEN * WIN_LEN) as f32;
        //dbg!("cov_mat: {:?}", self.cov_mat);
    }

    pub fn cal_inner_mat(&mut self, eps: f32) {

        let mut mid = FixMatrix33::new();
        let win_size = WIN_LEN * WIN_LEN;
        mid[(0, 0)] = eps / (win_size as f32);
        mid[(1, 1)] = eps / (win_size as f32);
        mid[(2, 2)] = eps / (win_size as f32);
        mid += self.cov_mat;

        //println!("into for");
        for i in 0..win_size {
            for j in 0..win_size {
                //println!("doing {}", i as f32 / win_size as f32);
                let g_i = self.get_global_elem_idx(i);
                let g_j = self.get_global_elem_idx(j);

                if g_i > g_j {
                    continue;
                }

                let mut v = (self.mat_ref.elem[g_i] - self.mean_v).dot((mid * (self.mat_ref.elem[g_j] - self.mean_v).transpose()).transpose());

                v = (1.0 + v) / (-(win_size as f32));

                if g_i == g_j {
                    v = v + 1.0;
                }

                self.inner_mat[(i, j)] = v;
            }
        }

        println!("{:?}", self.inner_mat);
    }

    pub fn reduce(&self, map: &mut ItemMap) {
        //just consider the upper triangle(g_i <= g_j)
        let win_size = WIN_LEN * WIN_LEN;
        for i in 0..win_size {
            for j in i..win_size {
                let g_i = self.get_global_elem_idx(i);
                let g_j = self.get_global_elem_idx(j);

                if g_i <= g_j {
                    let mut v = 0.0;
                    match map.get(&(g_i, g_j)) {
                        Some(x) => v = *x,
                        None => {},
                    };

                    map.insert((g_i, g_j), v);
                }
            }
        }
    }

    pub fn get_global_elem_idx(&self, local_idx: usize) -> usize {
        let r = local_idx / WIN_LEN;
        let c = local_idx % WIN_LEN;

        (self.start_idx.0 + r) * self.mat_ref.ncol() + self.start_idx.1 + c
    }

    pub fn get_global_idx(&self, local_idx: usize) -> (usize, usize) {
        let r = local_idx / WIN_LEN;
        let c = local_idx % WIN_LEN;

        (self.start_idx.0 + r, self.start_idx.1 + c)
    }

    pub fn get_local_idx(&self, global_idx: usize) -> usize {
        let r = global_idx / self.mat_ref.ncol() - self.start_idx.0;
        let c = global_idx % self.mat_ref.ncol() - self.start_idx.1;
        //(r, c)
        r * WIN_LEN + c
    }
}


pub struct WinRefMatrix<'a> {
    elem:   Vec<WindowRef<'a>>,
    nrow:   usize,
    ncol:   usize,
}

impl<'a> WinRefMatrix<'a> {
    pub fn from_vec(elem: Vec<WindowRef<'a>>, nrow: usize, ncol: usize) -> WinRefMatrix<'a> {
        WinRefMatrix {
            elem: elem,
            nrow: nrow,
            ncol: ncol,
        }
    }
    pub fn from_rgb_matrix(rgb_mat: &'a RgbMatrix, win_row_count: usize, win_col_count: usize) -> WinRefMatrix<'a> {
        let mut v = Vec::<WindowRef<'a>>::with_capacity(win_row_count * win_col_count);
        

        for i in 0..win_row_count {
            for j in 0..win_col_count {
                v.push(WindowRef::<'a>::new(rgb_mat, (i, j)));
            }
        }

        WinRefMatrix {
            elem: v,
            nrow: win_row_count,
            ncol: win_col_count,
        }
    }

    pub fn from_rgb_matrix_parallel(rgb_matrix: &'a RgbMatrix) -> WinRefMatrix<'a> {
        unimplemented!();
    }

    pub fn cal_all_distribution(&mut self) {
        self.elem.iter_mut().for_each(|x| x.cal_distribution());
    }

    pub fn cal_all_distribution_parallel(&mut self) {
        unimplemented!();
    }

    pub fn cal_all_inner_mat(&mut self, eps: f32) {
        self.elem.iter_mut().for_each(|x| x.cal_inner_mat(eps));
    }

    pub fn cal_all_reduce(&mut self, map: &mut ItemMap) {
        self.elem.iter().for_each(|x| x.reduce(map));
    }
}






impl RgbMatrix {
    pub fn from_raw_vec(raw: Vec<u8>, nrow: usize, ncol: usize) -> RgbMatrix {
        let elem = (0..raw.len()).step_by(3).map(|i| {
            let mut c = Rgb::new();
            c[0] = raw[i] as f32 / 255.0;
            c[1] = raw[i + 1] as f32 / 255.0;
            c[2] = raw[i + 2] as f32 / 255.0;
            c
        }).collect();
        RgbMatrix { 
            elem:   elem,
            ncol:   ncol,
            nrow:   nrow,
        }
    }
    pub fn into_raw_vec() -> (Vec<u8>, usize, usize) {
        unimplemented!()
    }
}


impl<T: std::clone::Clone> DMatrix<T> {
    pub fn new(nrow: usize, ncol: usize) -> DMatrix<T> {
        let mut v = Vec::<T>::with_capacity(nrow * ncol);
        unsafe {v.set_len(nrow * ncol);}
        DMatrix {
            elem:   v,
            nrow:   nrow,
            ncol:   ncol,
        }
    }

    pub fn from_vec(v: Vec<T>, nrow: usize, ncol: usize) -> DMatrix<T> {
        assert_eq!(v.len(), nrow * ncol);
        DMatrix {
            elem:   v,
            nrow:   nrow,
            ncol:   ncol,
        }
    }

    pub fn into_vec(self) -> Vec<T> {
        self.elem
    }
    pub fn nrow(&self) -> usize {
        self.nrow
    }
    pub fn ncol(&self) -> usize {
        self.ncol
    }
}

macro_rules! matrix_fix_mxn_ty {
    //($name: ident, $elem_ty: ty, $init: expr, $len: expr) => {
    ($name: ident, $elem_ty: ty, $init: expr, $nrow: expr, $ncol: expr) => {
        #[derive(Clone, Copy)]
        pub struct $name {
            elem: [$elem_ty; $nrow * $ncol],
        }
        impl $name {
            pub fn new() -> $name {
                $name { elem: [$init; $nrow * $ncol] }
            }
            pub fn from_raw(elem: [$elem_ty; $nrow * $ncol]) -> $name {
                $name { elem: elem }
            }

            pub fn nrow(&self) -> usize {
                $nrow
            }
            pub fn ncol(&self) -> usize {
                $ncol
            }
            pub fn into_raw(self) -> [$elem_ty; $nrow * $ncol] {
                self.elem
            }
            pub fn into_vec(self) -> Vec<$elem_ty> {
                self.elem.to_vec()
            }

            pub fn clear(&mut self) {
                self.elem.iter_mut().for_each(|x| *x = $init);
            }
        }

        impl std::fmt::Debug for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                let mut s = String::new();
                for i in 0..self.nrow() {
                    for j in 0..self.ncol() {
                        s = s + &format!(" {}", self.elem[i * self.ncol() + j]);
                    }
                    s = s + "\n";
                }
                write!(f, "Matrix \n{}", s)
            }
        }

    }
}

macro_rules! matrix_fix_nxn_ty {
    ($name: ident, $elem_ty: ty, $init: expr, $len: expr) => {
        matrix_fix_mxn_ty!{$name, $elem_ty, $init, $len, $len}
    }
}


macro_rules! matrix_math_ops_for_single_elem {
    ($name: ident, $other_elem_ty: ty) => {
        impl std::ops::Add<$other_elem_ty> for $name {
            type Output = $name;
            fn add(self, rhs: $other_elem_ty) -> Self::Output {
                let mut m = self.clone();
                m.elem.iter_mut().for_each(|x| *x = *x + rhs);
                m
            }
        }

        impl std::ops::Sub<$other_elem_ty> for $name {
            type Output = $name;
            fn sub(self, rhs: $other_elem_ty) -> Self::Output {
                let mut m = self.clone();
                m.elem.iter_mut().for_each(|x| *x = *x - rhs);
                m
            }
        }

        impl std::ops::Mul<$other_elem_ty> for $name {
            type Output = $name;
            fn mul(self, rhs: $other_elem_ty) -> Self::Output {
                let mut m = self.clone();
                m.elem.iter_mut().for_each(|x| *x = *x * rhs);
                m
            }
        }

        impl std::ops::Div<$other_elem_ty> for $name {
            type Output = $name;
            fn div(self, rhs: $other_elem_ty) -> Self::Output {
                let mut m = self.clone();
                m.elem.iter_mut().for_each(|x| *x = *x / rhs);
                m
            }
        }

        impl std::ops::AddAssign<$other_elem_ty> for $name {
            fn add_assign(&mut self, other: $other_elem_ty) {
                self.elem.iter_mut().for_each(|x| *x = *x + other);
            }
        }        

        impl std::ops::SubAssign<$other_elem_ty> for $name {
            fn sub_assign(&mut self, other: $other_elem_ty) {
                self.elem.iter_mut().for_each(|x| * x = *x - other);
            }
        }

        impl std::ops::MulAssign<$other_elem_ty> for $name {
            fn mul_assign(&mut self, other: $other_elem_ty) {
                self.elem.iter_mut().for_each(|x| *x = *x * other);
            }
        }

        impl std::ops::DivAssign<$other_elem_ty> for $name {
            fn div_assign(&mut self, other: $other_elem_ty) {
                self.elem.iter_mut().for_each(|x| *x = *x / other);
            }
        }
    }
}


macro_rules! matrix_able_ops_for_whole {
    ($name: ident, $elem_ty: ty) => {
        impl std::ops::Add for $name {
            type Output = $name;
            fn add(self, rhs: $name) -> Self::Output {
                assert_eq!(self.nrow(), rhs.nrow());
                assert_eq!(self.ncol(), rhs.ncol());
                let mut m = self.clone();
                for ((elem, l), r) in m.elem.iter_mut().zip(&self.elem).zip(&rhs.elem) {
                    *elem = *l + *r;
                }
                m
            }
        }

        impl std::ops::Sub for $name {
            type Output = $name;
            fn sub(self, rhs: $name) -> Self::Output {
                assert_eq!(self.nrow(), rhs.nrow());
                assert_eq!(self.ncol(), rhs.ncol());
                let mut m = self.clone();
                for ((elem, l), r) in m.elem.iter_mut().zip(&self.elem).zip(&rhs.elem) {
                    *elem = *l - *r;
                }
                m
            }
        }

        impl std::ops::AddAssign for $name {
            fn add_assign(&mut self, other: Self) {
                self.elem.iter_mut().enumerate().for_each(|(i, x)| *x += other.elem[i]);
            }
        }        

        impl std::ops::SubAssign for $name {
            fn sub_assign(&mut self, other: Self) {
                self.elem.iter_mut().enumerate().for_each(|(i, x)| *x -= other.elem[i]);
            }
        }        
 
    }
}

macro_rules! matrix_product_op {
    ($name: ident, $rhs_ty: ty, $res_ty: ty) => {
        impl std::ops::Mul<$rhs_ty> for $name {
            type Output = $res_ty;
            fn mul(self, rhs: $rhs_ty) -> $res_ty {
                let mut m = <$res_ty>::new();
                assert_eq!(self.ncol(), rhs.nrow());
                assert_eq!(m.nrow(), self.nrow());
                assert_eq!(m.ncol(), rhs.ncol());
                for i in 0..m.nrow() {
                    for j in 0..m.ncol() {
                        for k in 0..self.ncol() {
                            m.elem[i * m.ncol() + j] += self.elem[i * self.ncol() + k] * rhs.elem[k * rhs.ncol() + j];
                        }
                    }
                }
                m
            }
        }
    }
}

macro_rules! matrix_double_idx {
    ($name: ident, $elem_ty: ty) => {
        impl std::ops::Index<(usize, usize)> for $name {
            type Output = $elem_ty;
            fn index(&self, idx: (usize, usize)) -> &Self::Output {
                &self.elem[idx.0 * self.ncol() + idx.1]
            }
        }

        impl std::ops::IndexMut<(usize, usize)> for $name {
            fn index_mut(&mut self, idx: (usize, usize)) -> &mut $elem_ty {
                let n = self.ncol();
                &mut self.elem[idx.0 * n + idx.1]
            }
        }
    }
}

macro_rules! matrix_single_idx {
    ($name: ident, $elem_ty: ty) => {
        impl std::ops::Index<usize> for $name {
            type Output = $elem_ty;
            fn index(&self, idx: usize) -> &Self::Output {
                &self.elem[idx]
            }
        }

        impl std::ops::IndexMut<usize> for $name {
            fn index_mut(&mut self, idx: usize) -> &mut $elem_ty {
                &mut self.elem[idx]
            }
        }
    }

}

macro_rules! matrix_mat_op {
    ($name: ident, $res_ty: ty) => {
        impl $name {
            //transpose
            pub fn transpose(self) -> $res_ty {
                let mut m = <$res_ty>::new();
                for i in 0..m.nrow() {
                    for j in 0..m.ncol() {
                        m.elem[i * m.ncol() + j] = self.elem[j * self.ncol() + i];
                    }
                }
                m
            }
        }
    }
}

macro_rules! matrix_reduce_op {
    ($name: ident, $elem_ty: ty, $init: expr) => {
        impl $name {
            pub fn sum(&self) -> $elem_ty {
                self.elem.iter().fold($init, |acc, x| acc + x)
            }
        }
    }
}


macro_rules! matrix_vec_dot_op {
    ($name: ident, $rhs_ty: ty, $res_ty: ty, $init: expr) => {
        impl $name {
            pub fn dot(&self, rhs: $rhs_ty) -> $res_ty {
                assert_eq!(self.elem.len(), rhs.elem.len());
                let mut s: $res_ty = $init;
                (0..self.elem.len()).for_each(|i| s += self.elem[i] * rhs.elem[i]);
                s
            }
        }
    }
}



//create_matrix
//matrix_fix_nxn_ty!{FixMatrix99, f32, 0.0, 9}
matrix_fix_nxn_ty!{FixMatrix33, f32, 0.0, 3}
matrix_fix_nxn_ty!{Window, f32, 0.0, WIN_LEN}
//create vector
matrix_fix_mxn_ty!{FixVec3, f32, 0.0, 1, 3}
//row vec 1 * 3
matrix_fix_mxn_ty!{FixRowVec3, f32, 0.0, 1, 3}
//col vec 3 * 1
matrix_fix_mxn_ty!{FixColVec3, f32, 0.0, 3, 1}




//math ops: add, sub. mul, div
matrix_math_ops_for_single_elem!{FixMatrix33, f32}
matrix_math_ops_for_single_elem!{Window, f32}
//matrix_math_ops_for_single_elem!{FixVec3, f32}
matrix_math_ops_for_single_elem!{FixRowVec3, f32}
matrix_math_ops_for_single_elem!{FixColVec3, f32}
matrix_math_ops_for_single_elem!{Matrix, f32}

//abel ops: add, sub
matrix_able_ops_for_whole!{RgbMatrix, Rgb}
//matrix_able_ops_for_whole!{FixVec3, f32}
matrix_able_ops_for_whole!{FixMatrix33, f32}
matrix_able_ops_for_whole!{FixRowVec3, f32}
matrix_able_ops_for_whole!{FixColVec3, f32}


//product
//3*3 \times 3*1
matrix_product_op!{FixMatrix33, FixColVec3, FixColVec3}
matrix_product_op!{FixColVec3, FixRowVec3, FixMatrix33}
//vec
matrix_vec_dot_op!{FixRowVec3, FixRowVec3, f32, 0.0}
matrix_vec_dot_op!{FixColVec3, FixColVec3, f32, 0.0}
//matrix_vec_dot_op!{FixRowVec3, FixColVec3, f32, 0.0}


//mat_op: transpose
matrix_mat_op!{FixRowVec3, FixColVec3}
matrix_mat_op!{FixColVec3, FixRowVec3}



//double idx api for matrix
matrix_double_idx!{FixMatrix33, f32}
matrix_double_idx!{Window, f32}
matrix_double_idx!{RgbMatrix, Rgb}
matrix_double_idx!{Matrix, f32}

//single idx api for vector
//matrix_single_idx!{FixVec3, f32}
matrix_single_idx!{FixRowVec3, f32}
matrix_single_idx!{FixColVec3, f32}


//reduce op: sum
matrix_reduce_op!{FixMatrix33, f32, 0.0}



impl FixMatrix33 {
    pub fn reverse(&self) -> FixMatrix33 {
        let d = self.det();
        let mut m = FixMatrix33::new();
        m[(0, 0)] = (self[(1, 1)] * self[(2, 2)] - self[(1, 2)] * self[(2, 1)]) / d;
        m[(0, 1)] = (self[(0, 2)] * self[(2, 1)] - self[(0, 1)] * self[(2, 2)]) / d;
        m[(0, 2)] = (self[(0, 1)] * self[(1, 2)] - self[(0, 2)] * self[(1, 1)]) / d;
        m[(1, 0)] = (self[(1, 2)] * self[(2, 0)] - self[(1, 0)] * self[(2, 2)]) / d;
        m[(1, 1)] = (self[(0, 0)] * self[(2, 2)] - self[(0, 2)] * self[(2, 0)]) / d;
        m[(1, 2)] = (self[(0, 2)] * self[(1, 0)] - self[(0, 0)] * self[(1, 2)]) / d;
        m[(2, 0)] = (self[(1, 0)] * self[(2, 1)] - self[(1, 1)] * self[(2, 0)]) / d;
        m[(2, 1)] = (self[(0, 1)] * self[(2, 0)] - self[(0, 0)] * self[(2, 1)]) / d;
        m[(2, 2)] = (self[(0, 0)] * self[(1, 1)] - self[(0, 1)] * self[(1, 0)]) / d;
        m
    }
    pub fn det(&self) -> f32 {
        self.elem[0] * (self.elem[4] * self.elem[8] - self.elem[5] * self.elem[7])
        - self.elem[3] * (self.elem[1] * self.elem[8] - self.elem[2] * self.elem[7])
        + self.elem[6] * (self.elem[1] * self.elem[5] - self.elem[2] * self.elem[1])
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reverse() {
        let mut m33 = FixMatrix33::new();
        m33[(0, 0)] = 5.0;
        m33[(1, 1)] = 2.0;
        m33[(2, 2)] = 1.0;

        let m = m33.reverse();
        println!("{:?}", m);

        println!("det = {}", m.det());

        println!("{:?}", m);

    }

}