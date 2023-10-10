use gb_io::seq::Feature;
use svg::node::element::{path::Data, Path, SVG};

use crate::dna_sequence::DNAsequence;

const FEATURE_HEIGHT: f64 = 100.0;
const FONT_SIZE_SMALL: i32 = 48;
const FONT_SIZE_LARGE: i32 = FONT_SIZE_SMALL*2;

#[derive(Clone, Debug)]
pub struct RenderSVG {
    pub document: svg::Document,
}

impl RenderSVG {
    pub fn from_dna_sequence(dna: &DNAsequence) -> Self {
        let c = svg::node::element::Circle::new()
            .set("cx",2000)
            .set("cy",1500)
            .set("r",1000)
            .set("fill","none")
            .set("stroke","black")
            .set("stroke-width",10)
            ;

        let label = svg::node::element::Text::new()
            .set("x",2000)
            .set("y",1500)
            .set("font-size",FONT_SIZE_LARGE)
            .set("font-family","Verdana")
            .set("text-anchor","middle")
            .set("fill","black")
            .add(svg::node::Text::new(dna.name().clone().unwrap_or_default()))
            ;

        let size_text = svg::node::element::Text::new()
            .set("x",2000)
            .set("y",1600)
            .set("font-size",FONT_SIZE_SMALL)
            .set("font-family","Verdana")
            .set("text-anchor","middle")
            .set("fill","grey")
            .add(svg::node::Text::new(format!("({} bp)",dna.len())))
            ;

        let mut ret = Self {
            document: svg::Document::new()
                .set("viewBox", (0, 0, 4000, 3000))
                .add(c)
                .add(label)
                .add(size_text)
        };

        let max_pos = dna.len() as i64;

        for feature in dna.features() {
            match feature.location {
                gb_io::seq::Location::Range(from, to) => {
                    ret.document = Self::render_feature_from_range(&mut ret, feature, max_pos, from,to);
                },
                gb_io::seq::Location::Between(_, _) => todo!(),
                gb_io::seq::Location::Complement(_) => todo!(),
                gb_io::seq::Location::Join(_) => todo!(),
                gb_io::seq::Location::Order(_) => todo!(),
                gb_io::seq::Location::Bond(_) => todo!(),
                gb_io::seq::Location::OneOf(_) => todo!(),
                gb_io::seq::Location::External(_, _) => todo!(),
                gb_io::seq::Location::Gap(_) => todo!(),
            }

        }

        ret
    }

    fn render_feature_from_range(ret: &RenderSVG, feature: &Feature, max_pos: i64, start: (i64,gb_io::seq::Before), end: (i64,gb_io::seq::After)) -> SVG {
        let from = start.0;
        let to = end.0;
        let inner = 1000.0 + (Self::feature_band(feature) as f64)*100.0 - FEATURE_HEIGHT/2.0 + FEATURE_HEIGHT/10.0;
        let outer = inner+FEATURE_HEIGHT*8.0/10.0;

        // Feature body
        let data = Data::new()
            .move_to(Self::pos2xy(from, max_pos, outer))
            .line_to(Self::pos2xy(from, max_pos, inner))
            .elliptical_arc_to(Self::pos2arc(to, max_pos, inner,1))
            .line_to(Self::pos2xy(to, max_pos, outer))
            .elliptical_arc_to(Self::pos2arc(from, max_pos, outer, 0))
            .close();
        let path = Path::new()
            .set("from",from)
            .set("to",to)
            .set("fill", Self::feature_color(feature))
            .set("stroke", "black")
            .set("stroke-width", 3)
            .set("d", data);

        let mut label_text = String::new();
        for k in ["gene","product","standard_name","protein_id"] {
            label_text = match feature.qualifier_values(k.into()).next() {
                Some(s) => s.to_owned(),
                None => continue,
            };
            break;
        }
        // eprintln!("{:?}",&feature.kv);

        // Feature label
        let label_color = match Self::feature_color(feature) {
            "white" => "black",
            other => other,
        };
        let band = Self::feature_band(feature) as f64;
        let (x,y) = Self::pos2xy((to+from)/2, max_pos, outer+FEATURE_HEIGHT*band);
        let label = svg::node::element::Text::new()
            .set("x",x)
            .set("y",y)
            .set("font-size",FONT_SIZE_SMALL)
            .set("font-family","Verdana")
            .set("text-anchor","middle")
            .set("fill",label_color)
            .add(svg::node::Text::new(label_text))
            ;

        ret.document.to_owned() // TODO FIXME
            .add(path)
            .add(label)
    }

    fn feature_band(feature: &Feature) -> i8 {
        match feature.kind.to_string().to_ascii_uppercase().as_str() {
            "CDS" => 1,
            "GENE" => -1,
            _ => 0
        }
    }

    fn feature_color(feature: &Feature) -> &str {
        match feature.kind.to_string().to_ascii_uppercase().as_str() {
            "CDS" => "red",
            "GENE" => "blue",
            _ => "white"
        }
    }

    fn pos2arc(pos: i64, max_pos: i64, radius: f64, flag: u8) -> (f64,f64,u8,u8,u8,f64,f64) {
        let (x,y) = Self::pos2xy(pos,max_pos,radius);
        (radius,radius,0,0,flag,x,y)
    }

    fn pos2xy(pos: i64, max_pos: i64, radius: f64) -> (f64,f64) {
        let degrees: f64 = 360.0 * (pos as f64) / (max_pos as f64) ;
        let t = degrees * std::f64::consts::PI / 180.0;
        let h = 2000.0;
        let k = 1500.0;
        let x = radius * t.cos() + h;
        let y = radius * t.sin() + k;
        (x,y)
    }
}