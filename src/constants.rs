use macroquad::prelude::*;

pub const DT: f32 = 1. / 60.;
pub const RADIUS: f32 = 3.;
pub const CELL_SIZE: f32 = 20.;

pub const OVER_RELAXATION: f32 = 1.9;
pub const DENSITY_STIFFNESS: f32 = 1.;
pub const FLIP_RATIO: f32 = 0.9;
pub const SEPARATION_STEPS: usize = 3;
pub const INCOMPLESSIBILITY_STEPS: usize = 5;
