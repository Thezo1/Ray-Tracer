#![allow(dead_code)]
use std::ops;
//(todo!) impl cross and dot multiplication
pub type Point = Vec3;
pub type Color = Vec3;
pub type Ray = Vec3;

#[derive(Debug, Clone, PartialEq, Copy)]
pub struct Vec3{
    pub e: [f32;3],
}

impl Vec3{
    pub fn new(e1: f32, e2: f32, e3: f32) -> Self{
        Self{
           e:  [e1, e2, e3],
        }
    }

    pub fn get(&self, coordinate: u8) -> f32{
        match coordinate{
            0 => self.e[0],
            1 => self.e[1],
            2 => self.e[2],
            _ => panic!("out of bounds"),
        }
    }

    pub fn mul(&self, multiplyer: f32) -> Self{
        Self::new(
            self.e[0] * multiplyer,
            self.e[1] * multiplyer,
            self.e[2] * multiplyer,
        )
    }

    
    pub fn div(&self, divisor: f32) -> Self{
        Self::new(
            self.e[0] / divisor,
            self.e[1] / divisor,
            self.e[2] / divisor,
        )
    }
    
    pub fn magnitude(&self) -> f32{
        //sqrt(self*self)
        let sum_of = (self.e[0]*self.e[0] )
        + (self.e[1]*self.e[1] )
        +(self.e[2]* self.e[2]);

        sum_of.sqrt()
    }

    pub fn to_unit_vector(self) -> Self{
        self.div(self.magnitude())
    }

    pub fn dot_multiplication(&self, other: Vec3) -> f32{
        let res = self.e[0] * other.e[0]+
                  self.e[1] * other.e[1]+
                  self.e[2] * other.e[2];
        return res;
    }

    pub fn cross_multiplication(&self, other: Vec3) -> Vec3{
        let new_x = self.e[1] * other.e[2] - self.e[2] * other.e[1];
        let new_y = self.e[0] * other.e[2] - self.e[2] * other.e[0];
        let new_z = self.e[0] * other.e[1] - self.e[1] * other.e[0];

        let res = Vec3::new(new_x, new_y, new_z);
        return res;
    }
}

impl ops::Add for Vec3{
    type Output = Self;
    fn add(self, other: Self) -> Self{
        Self{
            e: [(self.e[0] + other.e[0]), (self.e[1] + other.e[1]), (self.e[2] + other.e[2])] 
        }
    }
}

impl ops::Sub for Vec3{
    type Output = Self;
    fn sub(self, other: Self) -> Self{
        Self{
            e: [(self.e[0] - other.e[0]), (self.e[1] - other.e[1]), (self.e[2] - other.e[2])] 
        }
    }
}

impl ops::Mul for Vec3{
    type Output = Self;
    fn mul(self, other: Self) -> Self{
        Self{
            e: [(self.e[0] * other.e[0]), (self.e[1] * other.e[1]), (self.e[2] * other.e[2])] 
        }
    }
}

impl ops::Neg for Vec3{
    type Output = Self;
    fn neg(self) -> Self::Output{
        Self{
            e:[
                -(self.e[0]),
                -(self.e[1]),
                -(self.e[2]),
            ]
        }
    }
}


impl AsRef<Vec3> for Vec3{
    
    fn as_ref(&self) -> &Self{
        self
    }
}

impl AsMut<Vec3> for Vec3{
    
    fn as_mut(&mut self) -> &mut Self{
        self
    }
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn addition_works(){ 
        let  new_vec3 = Vec3::new(0.0, 0.0, 0.0);
        let  new_vec3_1 = Vec3::new(5.2, 2.0, 0.5);
        let  to_add = Vec3::new(3.0, 6.0, 8.0); 
        let  to_add_1 = Vec3::new(3.0, 6.0, 8.0); 

        let res_1 = new_vec3 + to_add;
        let res_2 = new_vec3_1 + to_add_1;

        assert_eq!(res_1.e, [3.0, 6.0, 8.0]); 
        assert_eq!(res_2.e, [8.2, 8.0, 8.5]); 
    }

    #[test]
    #[ignore]
    fn subtraction_works(){ 
        let  new_vec3 = Vec3::new(0.0, 0.0, 0.0);
        let  new_vec3_1 = Vec3::new(5.2, 2.0, 0.5);
        let  to_subtract = Vec3::new(3.0, 6.0, 8.0); 
        let  to_subtract_1 = Vec3::new(3.0, 6.0, 8.0); 

        let res_1 = new_vec3 - to_subtract;
        let res_2 = new_vec3_1 - to_subtract_1;

        assert_eq!(res_1.e, [-3.0, -6.0, -8.0]); 
        assert_eq!(res_2.e, [2.2, -4.0, -7.5]); 
    }

    fn multiplication_works(){ 
        let  new_vec3 = Vec3::new(0.0, 0.0, 0.0);
        let  new_vec3_1 = Vec3::new(5.2, 2.0, 0.5);
        let  to_multiply = Vec3::new(3.0, 6.0, 8.0); 
        let  to_multiply_1 = Vec3::new(3.0, 6.0, 8.0); 

        let res_1 = new_vec3 * to_multiply;
        let res_2 = new_vec3_1 * to_multiply_1;

        assert_eq!(res_1.e, [0.0, 0.0, 0.0]); 
        assert_eq!(res_2.e, [15.6, 12.0, 4.0]); 
    }

    #[test]
    fn neg_works(){ 
        let  new_vec3 = Vec3::new(5.2, 2.0, 0.5);
        let neg_vec3 = -new_vec3; 

        assert_eq!(neg_vec3.e[0], -5.2);
        assert_eq!(neg_vec3.e[1], -2.0);
        assert_eq!(neg_vec3.e[2], -0.5);
    }
    #[test]
    fn get_works(){ 
        let  new_vec3 = Vec3::new(5.2, 2.0, 0.5);
         
        let x = new_vec3.get(0) ;
        let y = new_vec3.get(1);
        let z = new_vec3.get(2);

        assert_eq!(x, 5.2); 
        assert_eq!(y, 2.0); 
        assert_eq!(z, 0.5); 
    }
    #[test]
    fn mul_works(){
        let  new_vec3 = Vec3::new(5.2, 2.0, 0.5);
        let res = new_vec3.mul(2.0);

        assert_eq!(res.e, [10.4, 4.0, 1.0]);
    }

    
    #[test]
    fn div_works(){
        let  new_vec3 = Vec3::new(5.2, 2.0, 0.5);
        let res = new_vec3.div(2.0);

        assert_eq!(res.e, [2.6, 1.0, 0.25]);
    }

    #[test]
    fn magnitude_works(){
        let  new_vec3 = Vec3::new(5.2, 2.0, 0.5);
        let res = new_vec3.magnitude() as i32;

        assert_eq!(res.abs(),5);
    }

    #[test]
    #[ignore]
    fn to_unit_works(){
        let  new_vec3 = Vec3::new(5.2, 2.0, 0.5);
        let changed = new_vec3.to_unit_vector();

        assert_eq!(changed.e, [0.93, 0.36, 0.09]);
    }

    #[test]
    #[ignore]
    fn cross_multiplication_works(){ 
        let  new_vec3 = Vec3::new(5.2, 2.0, 0.5);
        let  to_multiply = Vec3::new(3.0, 6.0, 8.0); 

        let res = new_vec3.cross_multiplication(to_multiply);

        assert_eq!(res.e, [13.0, 40.1, 25.2]); 
    }

}


