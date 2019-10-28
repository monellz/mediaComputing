#[macro_escape]

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


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn quadratic() {
        let a = [2.0, 4.0, 5.0];
        let b = [2.5, 1.0, 3.0];
        let M = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];

        let a = V3::from_raw(a);
        let b = V3::from_raw(b);
        let M = M33::from_raw(M);

        println!("{:?}",a);
        println!("{:?}",b);
        println!("{:?}",M);

        let rhs = M.dot(&b);
        let ans = [13.5, 33.0,52.5];
        let right_answer = V3::from_raw(ans);
        assert_eq!(right_answer, rhs);

        let res = M.quadratic(&a, &b);
        assert_eq!(res, 421.5);
    }
}