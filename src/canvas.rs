#![allow(dead_code)]
use crate::vec3::Color;

pub struct Canvas{
    pub width: u32,
    pub height: u32,
    pub background_color: Color,
    pub pixels: Vec<Color>,
}

impl Canvas{
    pub fn new(width: u32, height: u32, background_color: Color) -> Self{
        let mut vec: Vec<Color> = Vec::with_capacity((width * height) as usize);
        for _ in 0..width{
            for _ in 0..height{
                vec.push(background_color);
            }
        }
        return Self{
            width,
            height,
            background_color,
            pixels: vec,
        };
    }

    pub fn color_at(&mut self, x: i32, y: i32, color: Color){
        let screen_x = (self.width/2) as i32 + x;
        let screen_y = (self.height/2) as i32 - y;
        let index = calc_idx(screen_x, screen_y, self.width as i32);
        self.pixels[index] = color;
    }

    pub fn convert_to_ppm(&self){
        print!("P3\n{} {}\n255\n", self.width, self.height);
        for pixel in self.pixels.iter(){
            println!("{} {} {}", 
                pixel.e[0],
                pixel.e[1],
                pixel.e[2],
            );
        }
    }
}


pub fn calc_idx(x: i32, y: i32, width: i32) -> usize{
    let index = (y * width + x) as usize;
    return index;
}


#[cfg(test)]
pub mod tests{
    use super::*;
    
    #[test]
    fn create_new_canvas_works(){
        let color_green = Color::new(0.0, 255.0, 0.0);
        let canvas = Canvas::new(10, 20, color_green);

        assert_eq!(canvas.height, 20);
        assert_eq!(canvas.width, 10);
        assert_eq!(canvas.background_color, color_green);
        assert!(!canvas.pixels.is_empty());
    }

    #[test]
    fn color_at_works(){
        let color_green = Color::new(0.0, 255.0, 0.0);
        let color_red = Color::new(255.0, 0.0, 0.0);
        let mut canvas = Canvas::new(10, 20, color_green);
    
        canvas.color_at(0, 0, color_red);
        let idx = calc_idx(5, 10, canvas.width as i32);
        assert_eq!(canvas.pixels[idx], color_red);
    }
}

