use crate::vec3::*;
use crate::scene::*;

pub struct Ray{
    origin: Point,
    direction: Point,
}

impl Ray{
    pub fn new(origin: Point, direction: Point) -> Self{

        return Self{origin, direction};
    }

    pub fn position(&self, t: f32) -> Point{
        return self.origin + self.direction.mul(t);
    }

    pub fn intersect(&self, sphere: Sphere) -> Vec<f32>{
    }
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn build_ray_works(){
        let ray = Ray::new(Point::new(1., 2., 3.), Point::new(4., 5., 6.));

        let origin = Point::new(1., 2., 3.);
        let direction = Point::new(4., 5., 6.);

        assert_eq!(ray.origin, origin);
        assert_eq!(ray.direction, direction);
    }

    #[test]
    fn position_works(){
        let ray = Ray::new(Point::new(2., 3., 4.), Point::new(1., 0., 0.));
        let pos1 = ray.position(0.);
        let pos2 = ray.position(1.);
        let pos3 = ray.position(-1.);
        let pos4 = ray.position(2.5);

        assert_eq!(pos1, Point::new(2., 3., 4.));
        assert_eq!(pos2, Point::new(3., 3., 4.));
        assert_eq!(pos3, Point::new(1., 3., 4.));
        assert_eq!(pos4, Point::new(4.5, 3., 4.));
    }

    #[test]
    fn intersect_works(){
        let ray = Ray::new(Point::new(2., 3., 4.), Point::new(1., 0., 0.));
        let sphere = Sphere::new(50., Point::new(0.0, 0.0, 0.0), Color::new(0.0, 0.0, 0.0), 500);
        let intersect = ray.intersect(sphere);

        assert_eq!(intersect[0], 4.);
        assert_eq!(intersect[1], 6.);

    }
}
