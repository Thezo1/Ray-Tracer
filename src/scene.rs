#![allow(dead_code)]
use crate::vec3::*;
use crate::light::*;

#[derive(PartialEq, Clone, Debug)]
pub struct Sphere{
    pub r: f32,
    pub c: Point,
    pub color: Color,
    pub specular: i32,
}

pub struct Scene{
    pub spheres: Vec<Option<Sphere>>, 
    pub lights: Vec<Light>,
}

impl Sphere{
    pub fn new(r: f32, c: Point, color: Color, specular: i32) -> Self{
        Self{
            r,
            c,
            color,
            specular,
        }
    }
}

impl AsRef<Sphere> for Sphere{
    
    fn as_ref(&self) -> &Self{
        self
    }
}

impl Scene{
    pub fn new(spheres: Vec<Option<Sphere>>, lights: Vec<Light>) -> Self{
        Self{
            spheres,
            lights,
        }
    }
}
