use crate::constants::{CELL_SIZE, RADIUS};
use crate::grid::{CellType, MacGrid};
use crate::particle::Particle;
use macroquad::prelude::*;

pub fn render(particles: &[Particle], grid: &MacGrid) {
    for j in 0..grid.ny {
        for i in 0..grid.nx {
            if grid.cell_type[i + j * grid.nx] == CellType::Wall {
                let x = i as f32 * CELL_SIZE;
                let y = j as f32 * CELL_SIZE;
                let color = Color::new(0.5, 0.5, 0.5, 1.0);
                draw_rectangle(x, y, CELL_SIZE, CELL_SIZE, color);
            }
        }
    }

    for p in particles {
        draw_circle(p.pos.x, p.pos.y, RADIUS, SKYBLUE);
    }
}
