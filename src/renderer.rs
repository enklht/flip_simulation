use crate::constants::RADIUS;
use crate::particle::Particle;
use macroquad::prelude::*;

const BG_COLOR: Color = Color::new(0.2, 0.18, 0.18, 1.0);
const WATER_COLOR: Color = Color::new(0.3, 0.85, 0.95, 0.9);

pub fn render(particles: &[Particle]) {
    clear_background(BG_COLOR);

    for p in particles {
        draw_circle(p.pos.x, p.pos.y, RADIUS * 2., WATER_COLOR);
    }
}
