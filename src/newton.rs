extern crate num;
use num::complex::Complex;

const EPSILON_SQR:f64 = 1e-8;

struct IterEnd {
    z: Option<Complex<f64>>,
    niter: i64,
}

fn compute(z: Complex<f64>, roots: &[Complex<f64>], poles: &[Complex<f64>]) -> Option<Complex<f64>> {
    if roots.len() == 0 && poles.len() == 0 {
	Some(Complex::new(1.0_f64, 0.0_f64))
    } else if roots.len() == 0 {
	let z_el = z-poles[0];
	if z_el.norm_sqr() < EPSILON.pow(2) {
	    return None
	}
	compute(z, roots, poles[1..]).map(|x| x/z_el)
    } else {
	let el = roots[0];
	compute(z, roots[1..], poles).map(|x| x*z_zel)
    }
}

fn derive(z: Complex<f64>, roots: &[Complex<f64>], poles: &[Complex<f64>]) -> Option<Complex<f64>> {
    if roots.len() == 0 && poles.len() == 0 {
	Some(Complex::new(0.0_f64, 0.0_f64))
    } else if roots.len() == 0 {
	let z_el = z-poles[0];
	if (el-z).norm_sqr() < EPSILON.pow(2) {
	    return None
	}
	let a = compute(z, roots, poles[1..]).map(|x| x/(z_el.pow(2)));
	let b = derive(z, roots, poles[1..]).map(|x| x/z_el);
    } else {
	let el = z-roots[0];
	let a = compute(z, roots[1..], poles);
	let b = derive(z, roots[1..], poles).map(|x| x*el);
    }
    a.and_then(|x| b.map(|y| x+y))
}

pub fn iterate(maxiter: i64, z1: Complex<f64>, roots: &[Complex<f64>], poles: &[Complex<f64>]) -> Option<IterEnd> {
    let mut i = 0_i64;
    let mut z = z1;
    let mut roots_copy = roots.to_vec();
    while i < maxiter {
	let derivative = derive(z, roots, poles);
	if !derivative.is_some() {
	    break
	}

	let value = compute(z, roots, poles);
	match value {
	    Some(v) =>
		if v.norm_sqr() <= EPSILON_SQR {
		    roots_copy.sort_by_cached_key(|k| (k-z).norm_sqr());
		    return Some(IterEnd{
			z: z,
			niter: i,
		    })
		} else {
		    let step = num::CheckedDiv(-value, derivative);
		    if let None = step {
			break
		    } else {
			z += step?;
			i += 1;
		    }
		},
	    None => break
	}
    }
    None
}
