// basic svg drawing components

use crate::plot_svg::constants;
use std::default::Default;

// ----------------------------------------------------------------------------
// Text
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Font

pub struct Font {
    pub default: usvg::Font,
}

impl Default for Font {
    fn default() -> Self {
        Self::new()
    }
}

impl Font {
    pub fn new() -> Self {
        Font {
            default: usvg::Font {
                families: vec!["DejaVu Sans".to_string()],
                style: usvg::FontStyle::Normal,
                stretch: usvg::FontStretch::Normal,
                weight: 200,
            },
        }
    }
}

// ----------------------------------------------------------------------------
// TextDecoration

pub struct TextDecoration {
    pub default: usvg::TextDecoration,
}

impl Default for TextDecoration {
    fn default() -> Self {
        Self::new()
    }
}

impl TextDecoration {
    pub fn new() -> Self {
        TextDecoration {
            default: usvg::TextDecoration {
                underline: None,
                overline: None,
                line_through: None,
            },
        }
    }
}

// ----------------------------------------------------------------------------
// TextSpan

pub struct TextSpan {
    pub default: usvg::TextSpan,
}

impl Default for TextSpan {
    fn default() -> Self {
        Self::new()
    }
}

impl TextSpan {
    pub fn new() -> Self {
        TextSpan {
            default: usvg::TextSpan {
                start: 0,
                end: 1000,
                fill: Some(usvg::Fill {
                    paint: usvg::Paint::Color(usvg::Color::new_rgb(0, 0, 0)),
                    opacity: usvg::tiny_skia_path::NormalizedF32::new(1.).unwrap(),
                    rule: usvg::FillRule::NonZero,
                }),
                stroke: None,
                paint_order: usvg::PaintOrder::FillAndStroke,
                font: Font::default().default,
                font_size: usvg::NonZeroPositiveF32::new(constants::FONT_SIZE as f32)
                    .unwrap(),
                small_caps: false,
                apply_kerning: false,
                decoration: TextDecoration::default().default,
                dominant_baseline: usvg::DominantBaseline::Auto,
                alignment_baseline: usvg::AlignmentBaseline::Auto,
                baseline_shift: vec![usvg::BaselineShift::Baseline],
                visibility: usvg::Visibility::Visible,
                letter_spacing: 1.,
                word_spacing: 1.,
                text_length: Some(1.),
                length_adjust: usvg::LengthAdjust::Spacing,
            },
        }
    }
}

// ----------------------------------------------------------------------------
// CharacterPosition

pub struct CharacterPosition {
    pub default: usvg::CharacterPosition,
}

impl Default for CharacterPosition {
    fn default() -> Self {
        Self::new()
    }
}

impl CharacterPosition {
    pub fn new() -> Self {
        CharacterPosition {
            default: usvg::CharacterPosition {
                x: Some(0.),
                y: Some(0.),
                dx: Some(0.),
                dy: Some(0.),
            },
        }
    }
}

// ----------------------------------------------------------------------------
// Stroke
// ----------------------------------------------------------------------------

pub struct Stroke {
    pub default: usvg::Stroke,
}

impl Default for Stroke {
    fn default() -> Self {
        Self::new()
    }
}

impl Stroke {
    pub fn new() -> Self {
        Stroke {
            default: usvg::Stroke {
                paint: usvg::Paint::Color(usvg::Color::new_rgb(0, 0, 0)),
                dasharray: None,
                dashoffset: 0.,
                miterlimit: usvg::StrokeMiterlimit::new(1.),
                opacity: usvg::tiny_skia_path::NormalizedF32::new(1.).unwrap(),
                width: usvg::NonZeroPositiveF32::new(1.).unwrap(),
                linecap: usvg::LineCap::Butt,
                linejoin: usvg::LineJoin::Miter,
            },
        }
    }
}
