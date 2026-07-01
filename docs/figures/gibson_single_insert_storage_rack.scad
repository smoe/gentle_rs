// GENtle rack OpenSCAD export
// rack_id=rack-1
// template=storage_pcr_tube_rack
// family=storage
// format=pcr_tube_0_2ml
$fn = 56;

outer_width = 62.800;
outer_depth = 44.000;
rack_height = 12.000;
opening_diameter = 9.200;
floor_thickness = 1.200;
front_top_clearance = 3.000;
front_label_strip_depth = 6.000;
front_label_strip_recess = 0.800;
corner_radius = 1.400;

module gentle_rack() {
    difference() {
        cube([outer_width, outer_depth, rack_height], false);
        translate([6.400, 7.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([16.400, 7.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([26.400, 7.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([36.400, 7.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([46.400, 7.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([56.400, 7.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([6.400, 17.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([16.400, 17.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([26.400, 17.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([36.400, 17.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([46.400, 17.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([56.400, 17.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([6.400, 27.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([16.400, 27.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([26.400, 27.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([36.400, 27.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([46.400, 27.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([56.400, 27.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([6.400, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([16.400, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([26.400, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([36.400, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([46.400, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([56.400, 37.600, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([0.600, 0.800, rack_height - front_label_strip_recess]) cube([61.600, 1.400, front_label_strip_recess + 0.05], false);
    }
}

gentle_rack();
