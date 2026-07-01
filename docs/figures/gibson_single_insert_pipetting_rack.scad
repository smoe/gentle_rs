// GENtle rack OpenSCAD export
// rack_id=rack-1
// template=pipetting_pcr_tube_rack
// family=pipetting
// format=pcr_tube_0_2ml
$fn = 56;

outer_width = 87.200;
outer_depth = 60.200;
rack_height = 18.000;
opening_diameter = 9.200;
floor_thickness = 2.000;
front_top_clearance = 5.000;
front_label_strip_depth = 8.000;
front_label_strip_recess = 1.200;
corner_radius = 2.400;

module gentle_rack() {
    difference() {
        cube([outer_width, outer_depth, rack_height], false);
        translate([8.600, 9.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([22.600, 9.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([36.600, 9.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([50.600, 9.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([64.600, 9.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([78.600, 9.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([8.600, 23.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([22.600, 23.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([36.600, 23.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([50.600, 23.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([64.600, 23.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([78.600, 23.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([8.600, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([22.600, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([36.600, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([50.600, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([64.600, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([78.600, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([8.600, 51.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([22.600, 51.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([36.600, 51.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([50.600, 51.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([64.600, 51.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([78.600, 51.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([1.200, 0.800, rack_height - front_label_strip_recess]) cube([84.800, 3.400, front_label_strip_recess + 0.05], false);
    }
}

gentle_rack();
