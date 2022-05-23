use winter_math::StarkField;

/// Print field-element vector concisely as 0, 1, 2, * (for pos >2), - (for neg).
/// Example: ...1..1.1.2.*.-.1.. with no newline.
pub fn print_vec_bits<E: StarkField>(vec: &Vec<E>) {
    for elt in vec {
        if elt == &E::ZERO {
            print!(".");
        } else if elt == &E::ONE {
            print!("1");
        } else if elt == &E::ONE.neg() {
            print!("-");
        } else if elt == &E::from(2u64) {
            print!("2");
        } else {
            print!("*");
        }
    }
}

pub fn print_vec<E: StarkField>(vec: &Vec<E>) {
    for elt in vec {
        print!("{:?} ", elt.as_int());
    }
}
