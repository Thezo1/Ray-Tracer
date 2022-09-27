use crate::vec3::*;

#[derive(PartialEq, Clone)]
pub struct Light{
    pub ty: String,
    pub intensity: f32,
    pub stance: Option<Point>,
}

impl Light{
    pub fn new(ty: String, intensity: f32, stance: Option<Point>) ->Self{
        Self{
            ty,
            intensity,
            stance,
        }
    }
}

impl AsRef<Light> for Light{
    
    fn as_ref(&self) -> &Self{
        self
    }
}


