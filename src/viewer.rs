use crate::mesh_2d::Mesh2D;
use egui::{Color32, Pos2, Sense, Stroke};

/// Interactive 2D mesh viewer with optional scalar/vector overlays.
///
/// Display rules:
/// - If `triangle_field` is set, color each triangle by that value.
///   Otherwise color by the per-node average of `node_scalar`.
/// - If `displacement` is set, vertex world position = node + scale · disp.
///   The `show_undeformed` toggle additionally draws a faint wireframe at
///   the original (undeformed) positions.
///
/// Coordinate convention: world (mesh) coords use math-style Y up. Screen Y
/// goes down, so we flip Y when mapping world → screen.
pub struct FemViewer {
    pub mesh: Mesh2D,
    pub node_scalar: Option<Vec<f64>>,
    pub triangle_field: Option<Vec<f64>>,
    pub displacement: Option<Vec<[f64; 2]>>,
    pub displacement_scale: f64,

    pub umin: f64,
    pub umax: f64,
    pub show_wireframe: bool,
    pub show_undeformed: bool,

    center_x: f64,
    center_y: f64,
    scale: f32,
    fitted: bool,
}

impl FemViewer {
    pub fn new(mesh: Mesh2D) -> Self {
        Self {
            mesh,
            node_scalar: None,
            triangle_field: None,
            displacement: None,
            displacement_scale: 1.0,
            umin: 0.0,
            umax: 1.0,
            show_wireframe: true,
            show_undeformed: true,
            center_x: 0.0,
            center_y: 0.0,
            scale: 1.0,
            fitted: false,
        }
    }

    pub fn with_node_scalar(mut self, u: Vec<f64>) -> Self {
        if self.triangle_field.is_none() {
            self.umin = u.iter().copied().fold(f64::INFINITY, f64::min);
            self.umax = u.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        }
        self.node_scalar = Some(u);
        self
    }

    pub fn with_triangle_field(mut self, f: Vec<f64>) -> Self {
        self.umin = f.iter().copied().fold(f64::INFINITY, f64::min);
        self.umax = f.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        self.triangle_field = Some(f);
        self
    }

    pub fn with_displacement(mut self, d: Vec<[f64; 2]>, scale: f64) -> Self {
        self.displacement = Some(d);
        self.displacement_scale = scale;
        self
    }

    pub fn reset_view(&mut self) {
        self.fitted = false;
    }

    fn vertex_world(&self, node: usize) -> [f64; 2] {
        let p = self.mesh.nodes[node];
        match &self.displacement {
            Some(d) => [
                p[0] + self.displacement_scale * d[node][0],
                p[1] + self.displacement_scale * d[node][1],
            ],
            None => p,
        }
    }

    fn bbox(&self) -> (f64, f64, f64, f64) {
        (0..self.mesh.num_nodes()).fold(
            (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            |(xa, xb, ya, yb), n| {
                let [x, y] = self.vertex_world(n);
                (xa.min(x), xb.max(x), ya.min(y), yb.max(y))
            },
        )
    }

    fn fit_to_rect(&mut self, rect: egui::Rect) {
        let (xmin, xmax, ymin, ymax) = self.bbox();
        let bw = (xmax - xmin).max(1e-12) as f32;
        let bh = (ymax - ymin).max(1e-12) as f32;
        let sx = rect.width() / bw;
        let sy = rect.height() / bh;
        self.scale = sx.min(sy) * 0.9;
        self.center_x = 0.5 * (xmin + xmax);
        self.center_y = 0.5 * (ymin + ymax);
        self.fitted = true;
    }

    fn world_to_screen(&self, world: [f64; 2], canvas_center: Pos2) -> Pos2 {
        Pos2 {
            x: (world[0] - self.center_x) as f32 * self.scale + canvas_center.x,
            y: -(world[1] - self.center_y) as f32 * self.scale + canvas_center.y,
        }
    }

    fn screen_to_world(&self, screen: Pos2, canvas_center: Pos2) -> [f64; 2] {
        [
            ((screen.x - canvas_center.x) / self.scale) as f64 + self.center_x,
            (-(screen.y - canvas_center.y) / self.scale) as f64 + self.center_y,
        ]
    }

    pub fn show(&mut self, ui: &mut egui::Ui) {
        let avail = ui.available_size();
        let (response, painter) = ui.allocate_painter(avail, Sense::click_and_drag());
        let rect = response.rect;

        if !self.fitted {
            self.fit_to_rect(rect);
        }

        if response.dragged() {
            let d = response.drag_delta();
            self.center_x -= (d.x / self.scale) as f64;
            self.center_y += (d.y / self.scale) as f64;
        }

        if let Some(hover) = response.hover_pos() {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll.abs() > 0.0 {
                let canvas_center = rect.center();
                let before = self.screen_to_world(hover, canvas_center);
                let factor = (scroll * 0.0015).exp();
                self.scale = (self.scale * factor).clamp(1e-3, 1e6);
                let after = self.screen_to_world(hover, canvas_center);
                self.center_x += before[0] - after[0];
                self.center_y += before[1] - after[1];
            }
        }

        painter.rect_filled(rect, 0.0, Color32::from_gray(28));

        let canvas_center = rect.center();
        let wire = Stroke::new(0.5, Color32::from_black_alpha(80));
        let undeformed_stroke = Stroke::new(0.5, Color32::from_white_alpha(40));

        // Optionally draw the undeformed mesh as a faint wireframe behind.
        if self.show_undeformed && self.displacement.is_some() {
            for t in 0..self.mesh.num_triangles() {
                let [a, b, c] = self.mesh.triangles[t];
                let p0 = self.world_to_screen(self.mesh.nodes[a], canvas_center);
                let p1 = self.world_to_screen(self.mesh.nodes[b], canvas_center);
                let p2 = self.world_to_screen(self.mesh.nodes[c], canvas_center);
                painter.add(egui::Shape::line(
                    vec![p0, p1, p2, p0],
                    undeformed_stroke,
                ));
            }
        }

        for t in 0..self.mesh.num_triangles() {
            let [a, b, c] = self.mesh.triangles[t];
            let color_val = self.color_value(t, [a, b, c]);
            let color = self.colormap(color_val);
            let p0 = self.world_to_screen(self.vertex_world(a), canvas_center);
            let p1 = self.world_to_screen(self.vertex_world(b), canvas_center);
            let p2 = self.world_to_screen(self.vertex_world(c), canvas_center);
            painter.add(egui::Shape::convex_polygon(
                vec![p0, p1, p2],
                color,
                Stroke::NONE,
            ));
            if self.show_wireframe {
                painter.add(egui::Shape::line(vec![p0, p1, p2, p0], wire));
            }
        }
    }

    fn color_value(&self, t: usize, abc: [usize; 3]) -> f64 {
        if let Some(tf) = &self.triangle_field {
            tf[t]
        } else if let Some(u) = &self.node_scalar {
            let [a, b, c] = abc;
            (u[a] + u[b] + u[c]) / 3.0
        } else {
            0.0
        }
    }

    fn colormap(&self, v: f64) -> Color32 {
        let t = if (self.umax - self.umin).abs() < 1e-15 {
            0.5
        } else {
            ((v - self.umin) / (self.umax - self.umin)).clamp(0.0, 1.0)
        };
        if t < 0.5 {
            let s = (t * 2.0 * 255.0) as u8;
            Color32::from_rgb(s, s, 255)
        } else {
            let s = (255.0 - (t - 0.5) * 2.0 * 255.0) as u8;
            Color32::from_rgb(255, s, s)
        }
    }
}

/// An eframe::App wrapper around the viewer with a stats sidebar.
pub struct ViewerApp {
    pub label: String,
    pub viewer: FemViewer,
}

impl eframe::App for ViewerApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        egui::Panel::right("info").show_inside(ui, |ui| {
            ui.heading("FEM Viewer");
            ui.separator();
            ui.label(&self.label);
            ui.label(format!("nodes:     {}", self.viewer.mesh.num_nodes()));
            ui.label(format!("triangles: {}", self.viewer.mesh.num_triangles()));
            ui.label(format!(
                "boundary:  {}",
                self.viewer.mesh.boundary_nodes.len()
            ));
            ui.separator();
            let field_label = if self.viewer.triangle_field.is_some() {
                "stress / triangle field"
            } else {
                "u"
            };
            ui.label(format!("{field_label} min: {:+.4e}", self.viewer.umin));
            ui.label(format!("{field_label} max: {:+.4e}", self.viewer.umax));
            ui.separator();
            ui.checkbox(&mut self.viewer.show_wireframe, "show wireframe");
            if self.viewer.displacement.is_some() {
                ui.checkbox(&mut self.viewer.show_undeformed, "show undeformed");
                ui.horizontal(|ui| {
                    ui.label("disp scale:");
                    ui.add(
                        egui::Slider::new(
                            &mut self.viewer.displacement_scale,
                            0.0..=1000.0,
                        )
                        .logarithmic(true),
                    );
                });
            }
            if ui.button("reset view").clicked() {
                self.viewer.reset_view();
            }
            ui.separator();
            ui.label("drag to pan");
            ui.label("scroll to zoom");
        });

        egui::CentralPanel::default().show_inside(ui, |ui| {
            self.viewer.show(ui);
        });
    }
}
