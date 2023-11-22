use color_eyre::eyre::{Report, Result};

pub fn draw_raqote(
    canvas: &mut raqote::DrawTarget,
    x_coords: &[f32],
    y_coords: &[f32],
    fill: &raqote::Source,
    stroke: &raqote::Source,
    stroke_style: &raqote::StrokeStyle,
) -> Result<(), Report> {
    // construct a list of path operations (drawing instructions)
    let mut ops: Vec<raqote::PathOp> = Vec::new();
    let winding = raqote::Winding::EvenOdd;

    for (i, it) in x_coords.iter().zip(y_coords.iter()).enumerate() {
        let (x, y) = it;
        let point = raqote::Point::new(*x, *y);
        let op = match i {
            0 => raqote::PathOp::MoveTo(point),
            _ => raqote::PathOp::LineTo(point),
        };
        ops.push(op);
    }

    ops.push(raqote::PathOp::Close);

    let path = raqote::Path { ops, winding };

    // fill polygon
    canvas.fill(&path, fill, &raqote::DrawOptions::new());
    // outline polygon
    canvas.stroke(&path, stroke, stroke_style, &raqote::DrawOptions::new());

    Ok(())
}
