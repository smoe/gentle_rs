// GENtle rack OpenSCAD export
// rack_id=cell-culture-6well
// template=cell_culture_plate
// family=cell_culture
// format=cell_culture_plate_well
$fn = 56;

outer_width = 131.400;
outer_depth = 87.100;
rack_height = 4.000;
opening_diameter = 34.800;
floor_thickness = 1.000;
front_top_clearance = 4.000;
front_label_strip_depth = 8.000;
front_label_strip_recess = 0.800;
corner_radius = 2.000;

module gentle_rack() {
    difference() {
        cube([outer_width, outer_depth, rack_height], false);
        translate([26.600, 21.400, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([65.700, 21.400, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([104.800, 21.400, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([26.600, 60.500, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([65.700, 60.500, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([104.800, 60.500, -0.05]) cylinder(h=rack_height - floor_thickness + 0.10, d=opening_diameter);
        translate([0.500, 0.800, rack_height - front_label_strip_recess]) cube([130.400, 2.400, front_label_strip_recess + 0.05], false);
    }
}

gentle_rack();
