#ifndef PROBLEM_DERIVES_H
#define PROBLEM_DERIVES_H

// Step 3 evolution diagnostics -- see docs/teoria.md Sec 6.3/6.4 and
// Problem_Derive.H for what these mean and the arbitrary-z-axis caveat
// on the toroidal/poloidal split. E_pol is not stored separately: it
// is recovered in post-processing as emag_density - etor_density
// (fvolumesum on each, then subtract the volume integrals) rather than
// computed twice.
derive_lst.add("emag_density", IndexType::TheCellType(), 1, ca_deremag, the_same_box);
derive_lst.addComponent("emag_density", desc_lst, Mag_Type_x, 0, 1);
derive_lst.addComponent("emag_density", desc_lst, Mag_Type_y, 0, 1);
derive_lst.addComponent("emag_density", desc_lst, Mag_Type_z, 0, 1);

derive_lst.add("etor_density", IndexType::TheCellType(), 1, ca_deretor, the_same_box);
derive_lst.addComponent("etor_density", desc_lst, Mag_Type_x, 0, 1);
derive_lst.addComponent("etor_density", desc_lst, Mag_Type_y, 0, 1);

#endif
