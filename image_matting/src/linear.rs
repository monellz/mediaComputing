#[allow(dead_code)]
pub type Matrix = DMatrix<f64>;
#[derive(Debug, Clone)]
pub struct DMatrix<T> {
    pub elem: Vec<T>,
    nrow: usize,
    ncol: usize,
}


impl<T: std::clone::Clone> DMatrix<T> {
    pub fn new(nrow: usize, ncol: usize, init: T) -> DMatrix<T> {
        DMatrix {
            elem:   vec![init; nrow * ncol],
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

    pub fn copy_to_vec(&self) -> Vec<T> {
        self.elem.clone()
    }

    pub fn nrow(&self) -> usize {
        self.nrow
    }
    pub fn ncol(&self) -> usize {
        self.ncol
    }
}

#[derive(Clone, Copy)]
pub struct M33 {
    elem: [f64; 9],
}

impl M33 {
    pub fn eye(init: f64) -> M33 {
        let mut m = M33::new();
        m.elem[0] = init;
        m.elem[4] = init;
        m.elem[8] = init;
        m
    }
    pub fn from_v3_mul(lhs: &V3, rhs: &V3) -> M33 {
        let mut e = [0.0; 9];
        e[0] = lhs[0] * rhs[0];
        e[1] = lhs[0] * rhs[1];
        e[2] = lhs[0] * rhs[2];
        e[3] = lhs[1] * rhs[0];
        e[4] = lhs[1] * rhs[1];
        e[5] = lhs[1] * rhs[2];
        e[6] = lhs[2] * rhs[0];
        e[7] = lhs[2] * rhs[1];
        e[8] = lhs[2] * rhs[2];
        M33 { elem: e }
    }

    pub fn mat_mul(lhs: &M33, rhs: &M33) -> M33 {
        let mut e = [0.0; 9];
        e[0] = lhs.elem[0] * rhs.elem[0] + lhs.elem[1] * rhs.elem[3] + lhs.elem[2] * rhs.elem[6];
        e[1] = lhs.elem[0] * rhs.elem[1] + lhs.elem[1] * rhs.elem[4] + lhs.elem[2] * rhs.elem[7];
        e[2] = lhs.elem[0] * rhs.elem[2] + lhs.elem[1] * rhs.elem[5] + lhs.elem[2] * rhs.elem[8];
        e[3] = lhs.elem[3] * rhs.elem[0] + lhs.elem[4] * rhs.elem[3] + lhs.elem[5] * rhs.elem[6];
        e[4] = lhs.elem[3] * rhs.elem[1] + lhs.elem[4] * rhs.elem[4] + lhs.elem[5] * rhs.elem[7];
        e[5] = lhs.elem[3] * rhs.elem[2] + lhs.elem[4] * rhs.elem[5] + lhs.elem[5] * rhs.elem[8];
        e[6] = lhs.elem[6] * rhs.elem[0] + lhs.elem[7] * rhs.elem[3] + lhs.elem[8] * rhs.elem[6];
        e[7] = lhs.elem[6] * rhs.elem[1] + lhs.elem[7] * rhs.elem[4] + lhs.elem[8] * rhs.elem[7];
        e[8] = lhs.elem[6] * rhs.elem[2] + lhs.elem[7] * rhs.elem[5] + lhs.elem[8] * rhs.elem[8];

        M33 { elem: e }
    }

    pub fn nrow(&self) -> usize {
        3
    }
    pub fn ncol(&self) -> usize {
        3
    }

    pub fn quadratic(&self, lhs: &V3, rhs: &V3) -> f64 {
        //res = a Q b
        lhs[0] * (self.elem[0] * rhs[0] + self.elem[1] * rhs[1] + self.elem[2] * rhs[2])
        + lhs[1] * (self.elem[3] * rhs[0] + self.elem[4] * rhs[1] + self.elem[5] * rhs[2])
        + lhs[2] * (self.elem[6] * rhs[0] + self.elem[7] * rhs[1] + self.elem[8] * rhs[2])
    }

    pub fn reverse(&self) -> M33 {
        let d = self.det();
        let mut m = M33::new();
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

    pub fn det(&self) -> f64 {
        self.elem[0] * (self.elem[4] * self.elem[8] - self.elem[5] * self.elem[7])
        - self.elem[3] * (self.elem[1] * self.elem[8] - self.elem[2] * self.elem[7])
        + self.elem[6] * (self.elem[1] * self.elem[5] - self.elem[2] * self.elem[4])
    }

    pub fn dot(&self, rhs: &V3) -> V3 {
        let mut raw = [0.0; 3];
        raw[0] = self.elem[0] * rhs.elem[0] + self.elem[1] * rhs.elem[1] + self.elem[2] * rhs.elem[2];
        raw[1] = self.elem[3] * rhs.elem[0] + self.elem[4] * rhs.elem[1] + self.elem[5] * rhs.elem[2];
        raw[2] = self.elem[6] * rhs.elem[0] + self.elem[7] * rhs.elem[1] + self.elem[8] * rhs.elem[2];
        V3::from_raw(raw)
    }
}

impl std::fmt::Debug for M33 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "-----M33-----\n").unwrap();
        write!(f, "{}, {}, {}\n", self[(0, 0)], self[(0, 1)], self[(0, 2)]).unwrap();
        write!(f, "{}, {}, {}\n", self[(1, 0)], self[(1, 1)], self[(1, 2)]).unwrap();
        write!(f, "{}, {}, {}\n", self[(2, 0)], self[(2, 1)], self[(2, 2)]).unwrap();
        write!(f, "-------------\n")
    }
}



#[derive(Debug, Clone, Copy, PartialEq)]
pub struct V3 {
    elem: [f64; 3],
}

impl V3 {
    pub fn norm2(&self) -> f64 {
        self.elem[0] * self.elem[0] + self.elem[1] * self.elem[1] + self.elem[2] * self.elem[2]
    }
}

#[macro_export]
macro_rules! linear_basic_impl {
    ($name: ident, $elem_ty: ty, $init: expr, $elem_len: expr) => {
        impl $name {
            pub fn new() -> $name {
                $name { elem: [$init; $elem_len] }
            }
            pub fn from_raw(raw: [$elem_ty; $elem_len]) -> $name {
                $name { elem: raw }
            }
            pub fn to_raw(self) -> [$elem_ty; $elem_len] {
                self.elem
            }
            pub fn copy_to_raw(&self) -> [$elem_ty; $elem_len] {
                self.elem.clone()
            }

            pub fn copy_to_vec(&self) -> Vec<$elem_ty> {
                self.elem.to_vec()
            }
        }
    }
}


#[macro_export]
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
    };
}

#[macro_export]
macro_rules! matrix_able_ops_for_whole {
    ($name: ident, $elem_ty: ty) => {
        impl std::ops::Add for $name {
            type Output = $name;
            fn add(self, rhs: $name) -> Self::Output {
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

#[macro_export]
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

#[macro_export]
macro_rules! matrix_double_idx {
    ($name: ident, $elem_ty: ty, $nrow: expr, $ncol: expr) => {
        impl std::ops::Index<(usize, usize)> for $name {
            type Output = $elem_ty;
            fn index(&self, idx: (usize, usize)) -> &Self::Output {
                &self.elem[idx.0 * $ncol + idx.1]
            }
        }

        impl std::ops::IndexMut<(usize, usize)> for $name {
            fn index_mut(&mut self, idx: (usize, usize)) -> &mut $elem_ty {
                &mut self.elem[idx.0 * $ncol + idx.1]
            }
        }
    };
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

#[macro_export]
macro_rules! matrix_reduce_op {
    ($name: ident, $elem_ty: ty, $init: expr) => {
        impl $name {
            pub fn sum(&self) -> $elem_ty {
                self.elem.iter().fold($init, |acc, x| acc + x)
            }
        }
    }
}

linear_basic_impl!{M33, f64, 0.0, 9}
linear_basic_impl!{V3, f64, 0.0, 3}

matrix_math_ops_for_single_elem!{M33, f64}
matrix_math_ops_for_single_elem!{V3, f64}

matrix_able_ops_for_whole!{M33, f64}
matrix_able_ops_for_whole!{V3, f64}

matrix_single_idx!{V3, f64}
matrix_double_idx!{M33, f64, 3, 3}

matrix_reduce_op!{V3, f64, 0.0}


#[allow(dead_code)]
pub mod sparse {
    use num_traits;
    use rayon::prelude::*;

    pub type SpMatrix = CSRMatrix<f64>;
    pub type SpVec = SparseVector<f64>;

    pub struct SparseVector<T> {
        vals: Vec<T>,
        pos: Vec<usize>,
    }

    impl<T> SparseVector<T> 
        where T: std::marker::Copy + std::ops::AddAssign  {
        pub fn new(d: Vec<T>, p: Vec<usize>) -> SparseVector<T> {
            SparseVector {
                vals: d,
                pos: p,
            }
        }
    }


    pub struct CSRMatrix<T> {
        vals: Vec<T>,
        cols: Vec<usize>,
        row_offset: Vec<usize>,
    }
    impl<T> CSRMatrix<T>
        where T: std::marker::Copy + std::ops::AddAssign  {
        pub fn new(d: Vec<T>, r: Vec<usize>, c: Vec<usize>) -> CSRMatrix<T> {
            let mut rc: Vec<(usize, usize, T)> = Vec::with_capacity(d.len());
            assert_eq!(d.len(), r.len());
            assert_eq!(r.len(), c.len());

            for i in 0..d.len() {
                rc.push((r[i], c[i], d[i]));
            }
            rc.sort_by_key(|i| (i.0, i.1));

            let mut vals: Vec<T> = Vec::with_capacity(d.len() / 4 + 1);
            let mut cols: Vec<usize> = Vec::with_capacity(d.len() / 4 + 1);
            let mut row_offset: Vec<usize> = Vec::with_capacity(d.len() / 4 + 1);

            let mut row_cur_num = rc[0].0;            
            let mut offset_num = 0;
            let mut col_cur_num = rc[0].1;
            vals.push(rc[0].2);
            cols.push(rc[0].1);
            row_offset.push(0);

            for i in 1..rc.len() {
                if row_cur_num == rc[i].0 && col_cur_num == rc[i].1 {
                    vals[offset_num] += rc[i].2;
                } else {
                    if row_cur_num < rc[i].0 {
                        row_offset.push(offset_num + 1);
                    }
                    
                    vals.push(rc[i].2);
                    cols.push(rc[i].1);
                    
                    offset_num += 1;
                    row_cur_num = rc[i].0;
                    col_cur_num = rc[i].1;
                }
            }
            row_offset.push(offset_num + 1);

            CSRMatrix {
                vals,
                cols,
                row_offset,
            }
        }

        pub fn get_vals_ref(&self) -> &Vec<T> {
            &self.vals
        }

        pub fn get_row_offset_ref(&self) -> &Vec<usize> {
            &self.row_offset
        }

        pub fn get_cols_ref(&self) -> &Vec<usize> {
            &self.cols
        }
    }

    impl<T: std::fmt::Debug> std::fmt::Debug for CSRMatrix<T> {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "----CSR matrix----\n").unwrap();
            /*
            let mut rec = 0;
            let mat_len = self.row_offset.len() - 1;
            for r in 0..mat_len {
                for c in 0..mat_len {
                    if self.cols[rec] < c {
                        write!(f, "0 ").unwrap();
                    } else if self.cols[rec] == c {
                        write!(f, "{:?} ", self.vals[rec]).unwrap();
                        rec += 1;
                    }
                }
                write!(f, "\n").unwrap();
            }
            */
            println!("{:?}", self.vals);
            println!("{:?}", self.cols);
            println!("{:?}", self.row_offset);

            write!(f, "----CSR matrix---over---\n")
        }
    }

    impl<T> CSRMatrix<T>
        where T: num_traits::float::Float
        + std::fmt::Display
        + std::ops::AddAssign
        + std::iter::Sum
        + std::marker::Send
        + std::marker::Sync {
        pub fn mult_vec_parallel(&self, rhs: &Vec<T>) -> Vec<T> {
            assert_eq!(rhs.len(), self.row_offset.len() - 1);
            let mut res: Vec<T> = vec![T::zero(); rhs.len()];
            res.par_iter_mut().enumerate().for_each(|(i, x)| {
                for j in self.row_offset[i]..self.row_offset[i + 1] {
                    *x += self.vals[j] * rhs[self.cols[j]];
                }
            });
            res
        }

        pub fn solve_pcg_parallel(&self, rhs: &Vec<T>, tol: T, max_steps: usize, x: Option<Vec<T>>, cond_rev: CSRMatrix<T>) -> Vec<T> {
            //solve A x = b
            let mut x = x.unwrap_or(vec![T::zero(); rhs.len()]);
            let mut cur_step: usize = 0;
            
            let mut r = self.mult_vec_parallel(&x);
            r.par_iter_mut().enumerate().for_each(|(i, x)| *x = rhs[i] - *x);

            //let norm = |v: &Vec<T>| v.par_iter().cloned().reduce(|| T::zero(), |acc, x| acc + x * x);
            let dot = |a: &Vec<T>, b: &Vec<T>| a.par_iter().cloned().zip(b).fold(|| T::zero(), |acc, (x, &y)| acc + x * y).sum();

            let mut d = cond_rev.mult_vec_parallel(&r);
            let mut delta_new = dot(&r, &d);
            let delta_0 = delta_new;

            while cur_step < max_steps && delta_new > tol * tol * delta_0 {
                if cur_step % 100 == 0 {
                    info!("{}/{} step, err = {}", cur_step, max_steps, delta_new - tol * tol * delta_0);
                }
                let q = self.mult_vec_parallel(&d);
                let alpha = delta_new / dot(&d, &q);
                
                //x = x + alpha d
                x.par_iter_mut().zip(&d).for_each(|(i, j)| *i = *i + alpha * *j);

                if cur_step % 50 == 0 {
                    let t = self.mult_vec_parallel(&x);
                    r.par_iter_mut().enumerate().for_each(|(i, x)| *x = rhs[i] - t[i]);
                } else {
                    r.par_iter_mut().zip(&q).for_each(|(i, j)| *i = *i - alpha * *j);
                }

                let s = cond_rev.mult_vec_parallel(&r);
                let delta_old = delta_new;
                delta_new = dot(&r, &s);
                let beta = delta_new / delta_old;
                d.par_iter_mut().enumerate().for_each(|(i, v)| *v = s[i] + beta * *v);
                cur_step += 1;
            }
            info!("parallel pcg over, step: {}", cur_step);
            x
        }
    }

    impl<T> CSRMatrix<T>
        where T: std::marker::Copy + std::ops::AddAssign
        + num_traits::float::Float + num_traits::ops::inv::Inv<Output = T> {
        pub fn get_jacobi_cond_rev(&self) -> CSRMatrix<T> {
            let mat_len = self.row_offset.len() - 1;
            let mut data = vec![T::zero(); mat_len];
            let rows: Vec<usize> = (0..mat_len).collect();
            let cols: Vec<usize> = (0..mat_len).collect();
            for r in 0..mat_len {
                for c in self.row_offset[r]..self.row_offset[r + 1] {
                    if self.cols[c] == r {
                        data[r] = self.vals[c].inv();
                    }
                }
            }
            CSRMatrix::new(data, rows, cols)
        }
        pub fn get_euclid_norm_cond_rev(&self) -> CSRMatrix<T> {
            let mat_len = self.row_offset.len() - 1;
            let mut data = vec![T::zero(); mat_len];
            let rows = (0..mat_len).collect();
            let cols = (0..mat_len).collect();
            for r in 0..mat_len {
                for c in self.row_offset[r]..self.row_offset[r + 1] {
                    data[self.cols[c]] += self.vals[c] * self.vals[c];
                }
            }
            for i in 0..data.len() {
                data[i] = data[i].sqrt().inv();
            }
            CSRMatrix::new(data, rows, cols)
        }

        pub fn get_abs_norm_cond_rev(&self) -> CSRMatrix<T> {
            let mat_len = self.row_offset.len() - 1;
            let mut data = vec![T::zero(); mat_len];
            let rows = (0..mat_len).collect();
            let cols = (0..mat_len).collect();
            for r in 0..mat_len {
                for c in self.row_offset[r]..self.row_offset[r + 1] {
                    data[self.cols[c]] += self.vals[c].abs();
                }
            }
            for i in 0..data.len() {
                data[i] = data[i].inv();
            }
            CSRMatrix::new(data, rows, cols)
        }
    }


    impl<T> CSRMatrix<T>
        where T: std::marker::Copy 
        + core::default::Default
        + num_traits::float::Float
        + std::cmp::PartialOrd
        + std::fmt::Display
        + std::ops::AddAssign {
        pub fn mult_vec(&self, rhs: &Vec<T>) -> Vec<T> {
            assert_eq!(rhs.len(), self.row_offset.len() - 1);

            let mut res: Vec<T> = vec![T::default(); rhs.len()];
            let mat_len = self.row_offset.len() - 1;
            for i in 0..mat_len {
                for j in self.row_offset[i]..self.row_offset[i + 1] {
                    res[i] += self.vals[j] * rhs[self.cols[j]];
                } 
            }
            res
        }

        pub fn solve_pcg(&self, rhs: &Vec<T>, tol: T, max_steps: usize, x: Option<Vec<T>>, cond_rev: CSRMatrix<T>) -> Vec<T> {
             //solve A x = b
            let mut x = x.unwrap_or(vec![T::default(); rhs.len()]);
            let mut cur_step: usize = 0;
            
            let mut r = self.mult_vec(&x);
            r.iter_mut().enumerate().for_each(|(i, x)| *x = rhs[i] - *x);

            let dot = |a: &Vec<T>, b: &Vec<T>| a.iter().zip(b).fold(T::default(), |acc, (&x, &y)| acc + x * y);

            let mut d = cond_rev.mult_vec(&r);
            let mut delta_new = dot(&r, &d);
            let delta_0 = delta_new;

            while cur_step < max_steps && delta_new > tol * tol * delta_0 {
                if cur_step % 100 == 0 {
                    info!("{}/{} step, err = {}", cur_step, max_steps, delta_new - tol * tol * delta_0);
                }
                let q = self.mult_vec(&d);
                let alpha = delta_new / dot(&d, &q);
                
                //x = x + alpha d
                x.iter_mut().zip(&d).for_each(|(i, j)| *i = *i + alpha * *j);

                if cur_step % 50 == 0 {
                    let t = self.mult_vec(&x);
                    r.iter_mut().enumerate().for_each(|(i, x)| *x = rhs[i] - t[i]);
                } else {
                    r.iter_mut().zip(&q).for_each(|(i, j)| *i = *i - alpha * *j);
                }

                let s = cond_rev.mult_vec(&r);
                let delta_old = delta_new;
                delta_new = dot(&r, &s);
                let beta = delta_new / delta_old;
                d.iter_mut().enumerate().for_each(|(i, v)| *v = s[i] + beta * *v);
                cur_step += 1;
            }
            info!("pcg over, step: {}", cur_step);
            x
        }
        pub fn solve_cg(&self, rhs: &Vec<T>, tol: T, max_steps: usize, x: Option<Vec<T>>) -> Vec<T> {
            //solve A x = b
            let mut x = x.unwrap_or(vec![T::default(); rhs.len()]);
            let mut cur_step: usize = 0;
            
            let mut r = self.mult_vec(&x);
            r.iter_mut().enumerate().for_each(|(i, x)| *x = rhs[i] - *x);

            let norm = |v: &Vec<T>| v.iter().fold(T::default(), |acc, &x| acc + x * x);
            let dot = |a: &Vec<T>, b: &Vec<T>| a.iter().zip(b).fold(T::default(), |acc, (&x, &y)| acc + x * y);

            let mut d = r.clone();
            let mut delta_new = norm(&r);
            let delta_0 = delta_new;

            while cur_step < max_steps && delta_new > tol * tol * delta_0 {
                if cur_step % 100 == 0 {
                    info!("{}/{} step, err = {}", cur_step, max_steps, delta_new - tol * tol * delta_0);
                }
 
                let q = self.mult_vec(&d);
                let alpha = delta_new / dot(&d, &q);
                
                //x = x + alpha d
                x.iter_mut().zip(&d).for_each(|(i, j)| *i = *i + alpha * *j);

                if cur_step % 50 == 0 {
                    let t = self.mult_vec(&x);
                    r.iter_mut().enumerate().for_each(|(i, x)| *x = rhs[i] - t[i]);
                } else {
                    r.iter_mut().zip(&q).for_each(|(i, j)| *i = *i - alpha * *j);
                }
                
                let delta_old = delta_new;
                delta_new = norm(&r);
                let beta = delta_new / delta_old;
                d.iter_mut().enumerate().for_each(|(i, v)| *v = r[i] + beta * *v);
                cur_step += 1;
            }
            info!("cg over, step: {}", cur_step);
            x
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use super::sparse::*;
    #[test]
    fn quadratic() {
        let a = [2.0, 4.0, 5.0];
        let b = [2.5, 1.0, 3.0];
        let mat = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];

        let a = V3::from_raw(a);
        let b = V3::from_raw(b);
        let mat = M33::from_raw(mat);

        println!("{:?}",a);
        println!("{:?}",b);
        println!("{:?}",mat);

        let rhs = mat.dot(&b);
        let ans = [13.5, 33.0,52.5];
        let right_answer = V3::from_raw(ans);
        assert_eq!(right_answer, rhs);

        let res = mat.quadratic(&a, &b);
        assert_eq!(res, 421.5);
    }

    #[test]
    fn sparse_mat_print() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        let rows = vec![0, 1, 2, 3];
        let cols = vec![0, 1, 2, 3];

        let mat = CSRMatrix::new(data, rows, cols);
        assert_eq!(mat.get_vals_ref(), &[1.0, 2.0, 3.0, 4.0]);
        assert_eq!(mat.get_cols_ref(), &[0, 1, 2, 3]);
        assert_eq!(mat.get_row_offset_ref(), &[0, 1, 2, 3, 4]);
        println!("{:?}", mat);

        /*
        1   7   0   0
        0   2   8   0
        5   0   3   9 
        0   6   0   4 
        */
        let data = vec![1.0, 7.0, 2.0, 8.0, 5.0, 3.0, 9.0, 6.0, 4.0];
        let rows = vec![0, 0, 1, 1, 2, 2, 2, 3, 3];
        let cols = vec![0, 1, 1, 2, 0, 2, 3, 1, 3];
        let mat = CSRMatrix::new(data, rows, cols);
        assert_eq!(mat.get_vals_ref(), &[1.0, 7.0, 2.0, 8.0, 5.0, 3.0, 9.0, 6.0, 4.0]);
        assert_eq!(mat.get_cols_ref(), &[0, 1, 1, 2, 0, 2, 3, 1, 3]);
        assert_eq!(mat.get_row_offset_ref(), &[0, 2, 4, 7, 9]);
        println!("{:?}", mat);

        let data = vec![1.0, 7.0, 2.0, 4.0, 5.0, 3.0, 9.0, 6.0, 4.0, 4.0];
        let rows = vec![0, 0, 1, 1, 2, 2, 2, 3, 3, 1];
        let cols = vec![0, 1, 1, 2, 0, 2, 3, 1, 3, 2];
        let mat = CSRMatrix::new(data, rows, cols);
        assert_eq!(mat.get_vals_ref(), &[1.0, 7.0, 2.0, 8.0, 5.0, 3.0, 9.0, 6.0, 4.0]);
        assert_eq!(mat.get_cols_ref(), &[0, 1, 1, 2, 0, 2, 3, 1, 3]);
        assert_eq!(mat.get_row_offset_ref(), &[0, 2, 4, 7, 9]);
        println!("{:?}", mat);
    }

    #[test]
    fn sparse_mat_mult_vec() {
        let data = vec![1.0, 7.0, 2.0, 4.0, 5.0, 3.0, 9.0, 6.0, 4.0, 4.0];
        let rows = vec![0, 0, 1, 1, 2, 2, 2, 3, 3, 1];
        let cols = vec![0, 1, 1, 2, 0, 2, 3, 1, 3, 2];
        let mat = CSRMatrix::new(data, rows, cols);
        
        let rhs = vec![1.0, 2.0, 3.0, 4.0];

        let res = mat.mult_vec(&rhs);
        println!("mult vec = {:?}", res);
        assert_eq!(res, vec![15.0, 28.0, 50.0, 28.0]);
    }

    #[test]
    fn sparse_solve() {
        /*
        1    0    5    0
        0    2    8    0
        5    8    3    9
        0    0    9    4
        */
        let data = vec![1.0, 5.0, 2.0, 8.0, 5.0, 8.0, 3.0, 9.0, 9.0, 4.0];
        let rows = vec![0, 0, 1, 1, 2, 2, 2, 2, 3, 3];
        let cols = vec![0, 2, 1, 2, 0, 1, 2, 3, 2, 3];
        let mat = CSRMatrix::new(data, rows, cols);
 

        let rhs = vec![1.0, 2.0, 3.0, 4.0];

        let (res, cur_step) = mat.solve_cg(&rhs, 1e-10, 100, Some(vec![1.0; 4]));
        println!("{:?}", res);
        println!("{:?}", cur_step);
    }

    #[test]
    fn sparse_solve_diag() {
        /*
        1 0 0 0 
        0 2 0 0
        0 0 4 0 
        0 0 0 8
        */
        let data = vec![1.0, 2.0, 4.0, 8.0];
        let rows = vec![0, 1, 2, 3];
        let cols = vec![0, 1, 2, 3];
        let mat = CSRMatrix::new(data, rows, cols);

        let rhs = vec![1.0, 1.0, 1.0, 1.0];

        let (res, cur_step) = mat.solve_cg(&rhs, 1e-10, 100, None);
        println!("{:?}", res);
        println!("{:?}", cur_step);
    }
}