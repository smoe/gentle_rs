use svg::node::element::{path::Data, Path};

use crate::dna_sequence::{DNAsequence, Feature};

#[derive(Clone, Debug)]
pub struct RenderSVG {
    pub document: svg::Document,
}

impl RenderSVG {
    pub fn from_dna_sequence(dna: &DNAsequence) -> Self {
        let font_size_small = 48;
        let font_size_large = font_size_small*2;
        
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
            .set("font-size",font_size_large)
            .set("font-family","Verdana")
            .set("text-anchor","middle")
            .set("fill","black")
            .add(svg::node::Text::new(dna.name()))
            ;

        let mut ret = Self {
            document: svg::Document::new()
                .set("viewBox", (0, 0, 4000, 3000))
                .add(c)
                .add(label)
        };

        let feature_height = 100.0;
        let max_pos = dna.len() as i64;

        for feature in dna.features() {
            let from = feature.location.start;
            let to = feature.location.stop;
            let inner = 1000.0 + (Self::feature_band(feature) as f64)*100.0 - feature_height/2.0 + feature_height/10.0;
            let outer = inner+feature_height*8.0/10.0;

            // Feature body
            let data = Data::new()
                .move_to(Self::pos2xy(from, max_pos, outer))
                .line_to(Self::pos2xy(from, max_pos, inner))
                .elliptical_arc_to(Self::pos2arc(to, max_pos, inner,1))
                .line_to(Self::pos2xy(to, max_pos, outer))
                .elliptical_arc_to(Self::pos2arc(from, max_pos, outer, 0))
                .close();
            let path = Path::new()
                .set("from",feature.location.start)
                .set("to",feature.location.stop)
                .set("fill", Self::feature_color(feature))
                .set("stroke", "black")
                .set("stroke-width", 3)
                .set("d", data);

            let mut label_text = String::new();
            for k in ["gene","product","standard_name","protein_id"] {
                label_text = match feature.kv.get(k) {
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
            let (x,y) = Self::pos2xy((to+from)/2, max_pos, outer+feature_height);
            let label = svg::node::element::Text::new()
                .set("x",x)
                .set("y",y)
                .set("font-size",font_size_small)
                .set("font-family","Verdana")
                .set("text-anchor","middle")
                .set("fill",label_color)
                .add(svg::node::Text::new(label_text))
                ;
                

            ret.document = ret.document
                .add(path)
                .add(label);
        }

        ret
    }

    fn feature_band(feature: &Feature) -> i8 {
        match feature.kind {
            crate::dna_sequence::FeatureKind::CDS => 1,
            crate::dna_sequence::FeatureKind::Gene => -1,
            _ => 0
        }
    }

    fn feature_color(feature: &Feature) -> &str {
        match feature.kind {
            crate::dna_sequence::FeatureKind::CDS => "red",
            crate::dna_sequence::FeatureKind::Gene => "blue",
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