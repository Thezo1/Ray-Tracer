#![allow(dead_code)]
use std::ops;
use std::convert::From;
use crate::vec3::Point;

pub type MatrixData = Vec<Vec<f32>>;

#[derive(Debug, PartialEq)]
pub struct Matrix{
    pub data: MatrixData,
    pub row: usize,
    pub col: usize,
}

impl Matrix{
    pub fn new(data: MatrixData) -> Self{
        let row = data.len(); 
        let col = data[0].len();
        return Self{ data, row, col};
    }

    pub fn matrix1x1_1x1_mul(&self, other: Matrix) -> Matrix{
        let m_data = vec![vec![0.0;1];4];
        let mut m = Matrix::new(m_data); 
        for row in 0..self.row{
            for col in 0..other.col{
                m.data[row][col] = self.data[row][0] * other.data[0][col]+
                    self.data[row][0] * other.data[1][col]+
                    self.data[row][0] * other.data[2][col]+
                    self.data[row][0] * other.data[3][col];
            }
        }

        return m;
    }

    pub fn matrix4x4_1x1_mul(&self, other: Matrix) -> Matrix{
        let m_data = vec![vec![0.0;1];4];
        let mut m = Matrix::new(m_data); 
        for row in 0..self.row{
            for col in 0..other.col{
                m.data[row][col] = self.data[row][0] * other.data[0][col]+
                    self.data[row][1] * other.data[1][col]+
                    self.data[row][2] * other.data[2][col]+
                    self.data[row][3] * other.data[3][col];
            }
        }

        return m;
    }
    pub fn mul_by_identity(&self) -> Matrix{
        let m_data = vec![vec![0.0;self.col];self.row];
        let  mut m = Matrix::new(m_data);
        let identity_data = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],];
        let identity = Matrix::new(identity_data);

        for row in 0..self.row{
            for col in 0..identity.col{
                m.data[row][col] = self.data[row][0] * identity.data[0][col]+
                    self.data[row][1] * identity.data[1][col]+
                    self.data[row][2] * identity.data[2][col]+
                    self.data[row][3] * identity.data[3][col];
            }
        }
        return m;
    }

    pub fn transpose_matrix(&self) -> Matrix{
        let m_data = vec![vec![0.0;self.row];self.col];
        let  mut m = Matrix::new(m_data);

        let mut row_mat = 0;
        let mut col_mat  = 0;
        for row in self.data.clone(){
            let temp = row.clone();
            for i in temp{
                m.data[row_mat][col_mat] = i;
                row_mat += 1;
            }
            row_mat -= row_mat;
            col_mat += 1;
        }

        return m;
    }

    pub fn determinat_2x2(&self) -> f32{
        return self.data[0][0] * self.data[1][1] - self.data[0][1] * self.data[1][0];
    }

    pub fn create_submatrix(&self, trunc_row: usize, trunc_col: usize) -> Matrix{
        let m_data = self.data.clone();
        let mut m = Matrix::new(m_data);

        let mut x = 0;
        while x < self.row{
            if x == trunc_row{
                m.data.remove(x);
                m.row -= 1;
                break;
            }
            x += 1;
        }

        let mut x1 = 0;
        while x1 < (self.col - 1){
            m.data[x1].remove(trunc_col);
            x1 += 1;
        }
        m.col -= 1;

        return m;
    }

    pub fn calc_minor(&self, row: usize, col: usize) -> f32{
        let sub_matrix = self.create_submatrix(row, col);
        return sub_matrix.determinat_2x2();
    }

    pub fn get_cofactor(&self, row: usize, col: usize) -> f32{
        let minor = self.calc_minor(row, col);
        if self.row == 4{
            if col == 0 || col == 1 || col == 2 || col == 3{
                let sub_matrix2 = self.create_submatrix(row, col);
                let mut determinant = sub_matrix2.determinant();
                   if (row + col)%2 != 0{
                      determinant = -determinant; 
                   }
                    return determinant;
            }else{
                let sub_matrix2 = self.create_submatrix(row, col - 1);
                let mut determinant = sub_matrix2.determinant();
                   if (row + col)%2 != 0{
                      determinant = -determinant; 
                   }
                    return determinant;
            }
        }
        if (row + col)%2 == 0{
            return minor;
        }else{
            return -minor;
        }
    }

    pub fn determinant(&self) -> f32{
        let mut determinant = 0.0;
        for x in 0..self.row{ 
            determinant += self.get_cofactor(0, x) * self.data[0][x];
        }

        return determinant;
    }

    pub fn can_calc_inverse(&self) -> bool{
        let mut res = true;
        if self.determinant() == 0.0{
            res = false;
        }
        return res;
    }
    
    pub fn calc_inverse(&self) -> Matrix{
        let m_data = vec![vec![0.0;self.row];self.col];
        let  mut m = Matrix::new(m_data);

        if !self.can_calc_inverse(){
            panic!("Cannot calc inverse, determinant is 0");
        }

        for row in 0..self.row {
            for col in 0..self.col{
                let cofactor = self.get_cofactor(row, col);

                m.data[col][row] = cofactor/self.determinant();
            }
        }

        return m;
    }
}

//Matrix Transformations
impl Matrix{
    ///set w = 1 for point,
    ///set w = 0 for vector,
    ///mul by a vector changes nothing
    fn get_translation_matrix(point: Point, w: f32) -> Matrix{

        let identity_data = vec![
            vec![1.0, 0.0, 0.0, point.e[0]],
            vec![0.0, 1.0, 0.0, point.e[1]],
            vec![0.0, 0.0, 1.0, point.e[2]],
            vec![0.0, 0.0, 0.0, w],];

        let m = Matrix::new(identity_data);

        return m;
    }

    pub fn translate(point_to_translate: Point, point: Point, w: f32) -> Point{
        let translation_matrix = Matrix::get_translation_matrix(point, w);

        let new_position = translation_matrix.matrix4x4_1x1_mul(Matrix::from(point_to_translate));

        let position: Point = Matrix::into(new_position);
        return position;
    }

    pub fn inverse_translate(point_to_translate: Point, point: Point, w: f32) -> Point{
        let translation_matrix = Matrix::get_translation_matrix(point, w);
        let inverse_translation = translation_matrix.calc_inverse();

        let new_position = inverse_translation.matrix4x4_1x1_mul(Matrix::from(point_to_translate));

        let position: Point = Matrix::into(new_position);
        return position;
    }

    fn get_scale_matrix(point: Point, w: f32) -> Matrix{
        let identity_data = vec![
            vec![point.e[0], 0.0, 0.0, 0.0],
            vec![0.0, point.e[1], 0.0, 0.0],
            vec![0.0, 0.0, point.e[2], 0.0],
            vec![0.0, 0.0, 0.0, w],];

        let m = Matrix::new(identity_data);

        return m;
    }

    pub fn scale(point_to_scale: Point, point:Point, w: f32) -> Point{
        let scale_matrix = Matrix::get_scale_matrix(point, w);

        let new_position = scale_matrix.matrix4x4_1x1_mul(Matrix::from(point_to_scale));

        let position: Point = Matrix::into(new_position);
        return position;
    }

    pub fn inverse_scale(point_to_scale: Point, point:Point, w: f32) -> Point{
        let scale_matrix = Matrix::get_scale_matrix(point, w);
        let inverse_scale = scale_matrix.calc_inverse();

        let new_position = inverse_scale.matrix4x4_1x1_mul(Matrix::from(point_to_scale));

        let position: Point = Matrix::into(new_position);
        return position;
    }

    pub fn rotate_x(point_to_rotate: Point, r: f32) -> Point{
        let transform_data = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, r.cos(), -r.sin(), 0.0],
            vec![0.0, r.sin(), r.cos(), 0.0],
            vec![0.0, 0.0, 0.0, 1.0],];

        let transform_matrix = Matrix::new(transform_data);
        let new_position = transform_matrix.matrix4x4_1x1_mul(Matrix::from(point_to_rotate));

        let position = Point::from(new_position);
        return position;
    }

    pub fn rotate_y(point_to_rotate: Point, r: f32) -> Point{
        let transform_data = vec![
            vec![r.cos(), 0.0, r.sin(), 0.0],
            vec![0.0, 1.0, 0.0,  0.0],
            vec![-r.sin(), 0.0,  r.cos(),  0.0],
            vec![0.0, 0.0, 0.0, 1.0],];

        let transform_matrix = Matrix::new(transform_data);
        let new_position = transform_matrix.matrix4x4_1x1_mul(Matrix::from(point_to_rotate));

        let position = Matrix::into(new_position);
        return position;
    }

    pub fn rotate_z(point_to_rotate: Point, r: f32) -> Point{
        let transform_data = vec![
            vec![r.cos(), -r.sin(), 0.0, 0.0],
            vec![r.sin(), r.cos(), 0.0,  0.0],
            vec![0.0, 0.0, 1.0,  0.0],
            vec![0.0, 0.0, 0.0, 1.0],];

        let transform_matrix = Matrix::new(transform_data);
        let new_position = transform_matrix.matrix4x4_1x1_mul(Matrix::from(point_to_rotate));

        let position = Matrix::into(new_position);
        return position;
    }

    pub fn inverse_rotate_x(point_to_rotate: Point, r: f32) -> Point{
        let transform_data = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, r.cos(), -r.sin(), 0.0],
            vec![0.0, r.sin(), r.cos(), 0.0],
            vec![0.0, 0.0, 0.0, 1.0],];

        let transform_matrix = Matrix::new(transform_data);
        let inverse = transform_matrix.calc_inverse();
        let new_position = inverse.matrix4x4_1x1_mul(Matrix::from(point_to_rotate));

        let position = Point::from(new_position);
        return position;
    }

    pub fn inverse_rotate_y(point_to_rotate: Point, r: f32) -> Point{
        let transform_data = vec![
            vec![r.cos(), 0.0, r.sin(), 0.0],
            vec![0.0, 1.0, 0.0,  0.0],
            vec![-r.sin(), 0.0,  r.cos(),  0.0],
            vec![0.0, 0.0, 0.0, 1.0],];

        let transform_matrix = Matrix::new(transform_data);
        let inverse = transform_matrix.calc_inverse();
        let new_position = inverse.matrix4x4_1x1_mul(Matrix::from(point_to_rotate));

        let position = Point::from(new_position);
        return position;
    }

    pub fn inverse_rotate_z(point_to_rotate: Point, r: f32) -> Point{
        let transform_data = vec![
            vec![r.cos(), -r.sin(), 0.0, 0.0],
            vec![r.sin(), r.cos(), 0.0,  0.0],
            vec![0.0, 0.0, 1.0,  0.0],
            vec![0.0, 0.0, 0.0, 1.0],];

        let transform_matrix = Matrix::new(transform_data);
        let inverse = transform_matrix.calc_inverse();
        let new_position = inverse.matrix4x4_1x1_mul(Matrix::from(point_to_rotate));

        let position = Point::from(new_position);
        return position;
    }

    pub fn shearing(point: Point, xy: f32, xz: f32, yx: f32, yz: f32, zx: f32, zy: f32) -> Point{
        let transform_data = vec![
            vec![1.0, xy, xz, 0.0],
            vec![yx, 1.0, yz, 0.0],
            vec![zx, zy, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],];

        let transform_matrix = Matrix::new(transform_data);
        let new_position = transform_matrix.matrix4x4_1x1_mul(Matrix::from(point));

        let position = Point::from(new_position);

        return position;
    }
}

impl From<Point> for Matrix{
    fn from(point: Point) -> Self{
        let m_data = vec![
            vec![point.e[0]],
            vec![point.e[1]],
            vec![point.e[2]],
            vec![1.],
        ];
        let m = Matrix::new(m_data);

        return m;
    }
}

impl From<Matrix> for Point{
    fn from(matrix: Matrix) -> Self{
        let point = Point::new(
            matrix.data[0][0],
            matrix.data[1][0],
            matrix.data[2][0],
        );

        return point;
    }
}

impl ops::Mul for Matrix{
    type Output = Self;
    fn mul(self, other: Self) -> Self{
        let m_data: Vec<Vec<f32>>= vec![vec![0.0;self.col];self.row];
        let mut m = Matrix::new(m_data);

        for row in 0..self.row{
            for col in 0..other.col{
                m.data[row][col] = self.data[row][0] * other.data[0][col]+
                    self.data[row][1] * other.data[1][col]+
                    self.data[row][2] * other.data[2][col]+
                    self.data[row][3] * other.data[3][col];
            }
        }
        return m;
    }
}

#[cfg(test)]
mod tests{
    use super::*;

   #[test]
   fn build_matrix_works(){
       let data = vec![vec![-3., 1.], vec![5., -2.]];
       let matrix = Matrix::new(data);

       assert_eq!(matrix.data[0][0], -3.);
       assert_eq!(matrix.data[0][1], 1.);
       assert_eq!(matrix.data[1][0], 5.);
       assert_eq!(matrix.data[1][1], -2.);
   }

   #[test]
   fn mul_matrices_work(){
       let m1_data =vec![
           vec![1., 2., 3., 4.],
          vec![5., 6., 7., 8.],
          vec![9., 8., 7., 6.],
          vec![5., 4., 3., 2.],];

       let m2_data = vec![
          vec![-2., 1., 2., 3.],
          vec![3., 2., 1., -1.],
          vec![4., 3., 6., 5.],
          vec![1., 2., 7., 8.],
       ];

       let res_data = vec![
          vec![20., 22., 50., 48.],
          vec![44., 54., 114., 108.],
          vec![40., 58., 110., 102.],
          vec![16., 26., 46., 42.],
       ];

       let m1 = Matrix::new(m1_data);
       let m2 = Matrix::new(m2_data);
       let res = Matrix::new(res_data);

       assert_eq!((m1 * m2), res);
   }

   #[test]
   fn matrix4x4_1x1_mul_works(){
       let m1_data =vec![
           vec![1., 2., 3., 4.],
           vec![2., 4., 4., 2.],
           vec![8., 6., 4., 1.],
           vec![0., 0., 0., 1.],];

       let m2_data = vec![
           vec![1.],
           vec![2.],
           vec![3.],
           vec![1.],];


       let res_data = vec![
           vec![18.],
           vec![24.],
           vec![33.],
           vec![1.],
       ];

       let m1 = Matrix::new(m1_data);
       let m2 = Matrix::new(m2_data);
       let res = Matrix::new(res_data);

       assert_eq!(m1.matrix4x4_1x1_mul(m2), res); 
   }

   #[test]
   fn mul_by_identity_works(){
       let m1_data =vec![
           vec![1., 2., 3., 4.],
           vec![2., 4., 4., 2.],
           vec![8., 6., 4., 1.],
           vec![0., 0., 0., 1.],];

       let m = Matrix::new(m1_data);

       let mul = m.mul_by_identity();

       assert_eq!(m, mul);

   }

   #[test]
   fn transpose_matrix_works(){
       let m1_data =vec![
           vec![1., 2., 3., 4.],
           vec![2., 4., 4., 2.],
           vec![8., 6., 4., 1.],
           vec![0., 0., 0., 1.],];

       let m1_transposed =vec![
           vec![1., 2., 8., 0.],
           vec![2., 4., 6., 0.],
           vec![3., 4., 4., 0.],
           vec![4., 2., 1., 1.],];

       let m = Matrix::new(m1_data);
       let m_trans = Matrix::new(m1_transposed);

       let m_transposed = m.transpose_matrix();

       assert_eq!(m_transposed, m_trans);

   }

   #[test]
   #[ignore]
   fn deterinant_2x2_works(){
       let m_data = vec![vec![1., 5.], vec![-3., 2.]];
       let m = Matrix::new(m_data);
       let deter = m.determinat_2x2();

       assert_eq!(deter, 17.);
   }

   #[test]
   fn calc_minor_works(){
       let m_data = vec![
           vec![3., 5., 0.],
           vec![2., -1., -7.],
           vec![6., -1., 5.],];

       let m = Matrix::new(m_data);
       let det_minor = m.calc_minor(1, 0);

       assert_eq!(det_minor, 25.);
   }

   #[test]
   fn get_cofactor_works(){
       let m_data = vec![
           vec![3., 5., 0.],
           vec![2., -1., -7.],
           vec![6., -1., 5.],];

       let m = Matrix::new(m_data);
       let cofactor = m.get_cofactor(0, 0);
       let cofactor1 = m.get_cofactor(1, 0);

       assert_eq!(cofactor, -12.);
       assert_eq!(cofactor1, -25.);
   }

   #[test]
   fn determinant_works(){
       let m_data = vec![
           vec![1., 2., 6.],
           vec![-5., 8., -4.],
           vec![2., 6., 4.],];
       
       let m2_data =vec![
           vec![-2., -8., 3., 5.],
           vec![-3., 1., 7., 3.],
           vec![1., 2., -9., 6.],
           vec![-6., 7., 7., -9.],];

       let m3_data =vec![
           vec![-5., 2., 6., -8.],
           vec![1., -5., 1., 8.],
           vec![7., 7., -6., -7.],
           vec![1., -3., 7., 4.],];


       let m4_data =vec![
           vec![6., 4., 4., 4.],
           vec![5., 5., 7., 6.],
           vec![4., -9., 3., -7.],
           vec![9., 1., 7., -6.],];

       let m = Matrix::new(m_data);
       let m2 = Matrix::new(m2_data);
       let m3 = Matrix::new(m3_data);
       let m4 = Matrix::new(m4_data);

       assert_eq!(m.determinant(), -196.);
       assert_eq!(m2.determinant(), -4071.);
       assert_eq!(m3.determinant(), 532.);
       assert_eq!(m4.determinant(), -2120.);
   }

   #[test]
   fn can_calc_inverse_works(){
       let m_data =vec![
           vec![6., 4.0, 4.0, 4.0],
           vec![5., 5.0, 7.0, 6.0],
           vec![4., -9.0, 3.0, -7.0],
           vec![9., 1.0, 7.0, -6.0],];

       let m2_data =vec![
           vec![-4., 2., -2., -3.],
           vec![9., 6., 2., 6.], vec![0., -5., 1., -5.],
           vec![0., 0., 0., 0.],];

       let m = Matrix::new(m_data);
       let m2 = Matrix::new(m2_data);

       assert_eq!(m.can_calc_inverse(), true);
       assert_eq!(m2.can_calc_inverse(), false);
   }

    #[test]
    #[ignore]
    fn calc_inverse_works(){
        let m_data =vec![
            vec![-5.0, 2.0, 6.0, -8.0],
            vec![1.0, -5.0, 1.0, 8.0],
            vec![7.0, 7.0, -6.0, -7.0],
            vec![1.0, -3.0, 7.0, 4.0],];

        let inverse_data =vec![
            vec![0.21805, 0.45113, 0.24060, 0.04511],
            vec![0.80827, -1.45677, 0.44361, 0.52068],
            vec![0.07895, 0.22368, 0.05263, 0.19737],
            vec![0.52256, 0.81391, 0.30075, 0.30639],];

        let m = Matrix::new(m_data);
        let inverse = Matrix::new(inverse_data);


        assert_eq!(m.determinant(), 532.0);
        assert_eq!(m.get_cofactor(2, 3), -160.0);
        assert_eq!(m.calc_inverse(), inverse);

    }
}

#[cfg(test)]
mod transformation_tests{
    use super::*;
    use std::f32::consts::PI;

   #[test]
   fn translation_works(){
       let translation_point = Point::new(5., -3., 2.);
       let to_be_translated = Point::new(-3., 4., 5.);
       let translate = Matrix::translate(to_be_translated, translation_point, 1.);

       let expected = Point::new(2., 1., 7.);

       assert_eq!(translate, expected);
   }

   #[test]
   fn inverse_translate_works(){
       let translation_point = Point::new(5., -3., 2.);
       let to_be_translated = Point::new(-3., 4., 5.);
       let translate = Matrix::inverse_translate(to_be_translated, translation_point, 1.);

       let expected = Point::new(-8., 7., 3.);

       assert_eq!(translate, expected);
   }

   #[test]
   fn scale_works(){
       let scale_point = Point::new(2., 3., 4.);
       let to_be_scaled = Point::new(-4., 6., 8.);
       let scale = Matrix::scale(to_be_scaled, scale_point, 1.);

       let vector = Point::new(-4., 6., 8.);
       let scale_vec = Matrix::scale(vector, scale_point, 0.);

       let expected = Point::new(-8., 18., 32.);

       assert_eq!(scale, expected);
       assert_eq!(scale_vec, expected);
   }
   //
   #[test]
   fn inverse_scale_works(){
       let scale_point = Point::new(2., 3., 4.);
       let to_be_scaled = Point::new(-4., 6., 8.);
       let scale = Matrix::inverse_scale(to_be_scaled, scale_point, 1.);

       let expected = Point::new(-2., 2., 2.);

       assert_eq!(scale, expected);
   }

   #[test]
   fn reflection_works(){
       let scale_point = Point::new(-1., 1., 1.);
       let to_be_scaled = Point::new(2., 3., 4.);
       let scale = Matrix::scale(to_be_scaled, scale_point, 1.);

       let expected = Point::new(-2., 3., 4.);

       assert_eq!(scale, expected);
   }

    #[test]
    #[ignore]
    fn rotation_x_works(){
        let point_to_rotate = Point::new(0., 1., 0.);
        let half_rotate = Matrix::rotate_x(point_to_rotate, PI/4.);
        let full_rotate = Matrix::rotate_x(point_to_rotate, PI/2.);

        let expected_half = Point::new(0., 2.2_f32.sqrt()/2., 2.2_f32.sqrt()/2.);
        let expected_full = Point::new(0., 0., 1.);

        assert_eq!(half_rotate, expected_half); 
        assert_eq!(full_rotate, expected_full); 
    }

    #[test]
    #[ignore]
    fn rotation_y_works(){
        let point_to_rotate = Point::new(0., 0., 1.);
        let full_rotate = Matrix::rotate_y(point_to_rotate, PI/2.);

        let expected_full = Point::new(1., 0., 0.);

        assert_eq!(full_rotate, expected_full); 
    }
    #[test]
    #[ignore]
    fn rotation_z_works(){
        let point_to_rotate = Point::new(0., 1., 0.);
        let full_rotate = Matrix::rotate_z(point_to_rotate, PI/2.);

        let expected_full = Point::new(-1., 0., 0.);

        assert_eq!(full_rotate, expected_full); 

    }

    #[test]
    #[ignore]
    fn inverse_rotation_x_works(){
        let point_to_rotate = Point::new(0., 1., 0.);
        let half_rotate = Matrix::inverse_rotate_x(point_to_rotate, PI/4.);

        let expected_half = Point::new(0., 2.2_f32.sqrt()/2., -2.2_f32.sqrt()/2.);

        assert_eq!(half_rotate, expected_half); 
    }

    #[test]
    fn shearing_works(){
        let point = Point::new(2., 3., 4.);

        let shearxy = Matrix::shearing(point, 1., 0., 0., 0., 0., 0.);
        let shearxz = Matrix::shearing(point, 0., 1., 0., 0., 0., 0.);
        let shearyx = Matrix::shearing(point, 0., 0., 1., 0., 0., 0.);
        let shearyz = Matrix::shearing(point, 0., 0., 0., 1., 0., 0.);
        let shearzx = Matrix::shearing(point, 0., 0., 0., 0., 1., 0.);
        let shearzy = Matrix::shearing(point, 0., 0., 0., 0., 0., 1.);

        let expected_shearxy = Point::new(5., 3., 4.);
        let expected_shearxz = Point::new(6., 3., 4.);
        let expected_shearyx = Point::new(2., 5., 4.);
        let expected_shearyz = Point::new(2., 7., 4.);
        let expected_shearzx = Point::new(2., 3., 6.);
        let expected_shearzy = Point::new(2., 3., 7.);

        assert_eq!(shearxy, expected_shearxy);
        assert_eq!(shearxz, expected_shearxz);
        assert_eq!(shearyx, expected_shearyx);
        assert_eq!(shearyz, expected_shearyz);
        assert_eq!(shearzx, expected_shearzx);
        assert_eq!(shearzy, expected_shearzy);
    }

    #[test]
    fn chaining_transformations(){
        let point = Point::new(1., 0., 1.);
        let expected = Point::new(15., 0., 7.);
        let scale_by =Point::new(5., 5., 5.);
        let translate_by =Point::new(10., 5., 7.);

        let rotate = Matrix::rotate_x(point, PI/2.);
        let scale = Matrix::scale(rotate, scale_by, 1.);
        let translate = Matrix::translate(scale, translate_by, 1.);

        assert_eq!(translate, expected);
    }
}

