#![allow(dead_code)]
use crate::vec3::Point;

pub type Vector = Point;

impl Vector{
    fn calc_start_point(&self, end_point: Point) -> Point{
        return end_point - *self;
    }

    fn calc_end_point(&self, start_point: Point) -> Point{
        return start_point + *self; 
    }

}


#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn calc_start_point_works(){
        let vec = Vector::new(2., 1., 0.);
        let end = Vector::new(5., 6., 0.);
        let start = Vector::new(3., 5., 0.);

        assert_eq!(vec.calc_start_point(end), start);
    }

    #[test]
    fn calc_end_point_works(){
        let vec = Vector::new(2., 1., 0.);
        let end = Vector::new(5., 6., 0.);
        let start = Vector::new(3., 5., 0.);

        assert_eq!(vec.calc_end_point(start), end);
    }
}
