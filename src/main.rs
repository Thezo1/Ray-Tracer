use std::f32::consts::PI;

use crate::vec3::*;
use crate::scene::*;
use crate::light::*;
use crate::vector::*;
use crate::ray::*;
use crate::canvas::*;
use crate::matrix::*;

mod vec3;
mod scene;
mod light;
mod vector;
mod canvas;
mod matrix; 
mod ray;


//Image Dimensions
const ASPECT_RATIO: f32 = 1.0/1.0;
const CAMERA_WIDTH: u32 = 400;
const CAMERA_HEIGHT: u32 = (CAMERA_WIDTH as f32 / ASPECT_RATIO) as u32;

//Camera Speccs
const VIEWPORT_HEIGHT: f32 = 1.0;
const VIEWPORT_WIDTH: f32 = ASPECT_RATIO * VIEWPORT_HEIGHT;
const FOCAL_LENGHT: f32 = 1.0;

fn main(){
    let  color_red: Color = Color::new(255.0, 0.0, 0.0);
    let  color_green: Color = Color::new(0.0, 255.0, 0.0);
    let  color_blue: Color = Color::new(0.0, 0.0, 255.0);
    let  color_white: Color = Color::new(255.0, 255.0, 255.0);
    let  color_black: Color = Color::new(0.0, 0.0, 0.0);

    let  mut canvas = Canvas::new(CAMERA_WIDTH, CAMERA_HEIGHT, color_blue);

    let center = Point::new(0., 0., 0.);
    let radius = 80.;

    let upper = center.e[1] + radius;

    let upper_point = Point::new(0., upper, 0.);

    let mut i=0;
    let mut pos = upper_point;
    while i < 12{
        let o_clock = Matrix::inverse_rotate_z(pos, PI/6.);
        canvas.color_at(o_clock.e[0] as i32, o_clock.e[1] as i32, color_white);
        pos = o_clock;
        i+=1;
    }

    
    canvas.convert_to_ppm();
}

