#[cfg(test)]
mod tests {
    use sprs;
    #[test]
    fn solver() {
        //3d
        let data = vec![2.0, 4.0, 8.0];
        let rows = vec![0, 1, 2];
        let cols = vec![0, 1, 2];

        let tri_mat = sprs::TriMat::from_triplets((3, 3), rows, cols, data);

        let tri_mat = tri_mat.to_csc();

        let mut rhs = vec![1.0, 1.0, 1.0];
        let res = sprs::linalg::trisolve::lsolve_csc_dense_rhs(tri_mat.view(), &mut rhs);

        assert_eq!(rhs, vec![0.5, 0.25, 0.125]);

        let copy = rhs.clone();
        assert_eq!(copy, vec![0.5, 0.25, 0.125]);
    }
}