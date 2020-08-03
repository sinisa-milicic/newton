extern crate num;
extern crate pyo3;


use num::complex::Complex;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use pyo3::types::*;
use pyo3::conversion::IntoPy;

type NewtComplex = Complex<f64>;
type NewtonResult = Result<NewtComplex, &'static str>;

const EPSILON_SQR:f64 = 1e-150;

pub struct IterEnd {
    z: Option<NewtComplex>,
    niter: i64,
}

impl IntoPy<PyObject> for IterEnd {
    fn into_py(self, _py: Python) -> PyObject {
	let mydict = PyDict::new(_py);
	mydict.set_item("z", self.z).unwrap();
	mydict.set_item("niter", self.niter).unwrap();
	(*mydict).to_object(_py)
    }
}

fn compute(z: NewtComplex, roots: &[NewtComplex],
	   poles: &[NewtComplex]) -> NewtonResult {
    if roots.len() == 0 && poles.len() == 0 {
	return Ok(Complex::new(1.0_f64, 0.0_f64))
    } else if roots.len() == 0 {
	let z_el = z-poles[0];
	if z_el.norm_sqr() < EPSILON_SQR.powf(2.0) {
	    return Err("Too close to pole")
	}
	compute(z, roots, &poles[1..]).map(|x| x/z_el)
    } else {
	let z_el = z-roots[0];
	compute(z, &roots[1..], poles).map(|x| x*z_el)
    }
}

fn derive(z: NewtComplex, roots: &[NewtComplex],
	  poles: &[NewtComplex]) -> NewtonResult {
    let a: NewtonResult;
    let b: NewtonResult;
    if roots.len() == 0 && poles.len() == 0 {
	return Ok(Complex::new(0.0_f64, 0.0_f64))
    } else if roots.len() == 0 {
	let z_el = z-poles[0];
	if (z_el).norm_sqr() < EPSILON_SQR.powf(2.0) {
	    return Err("Too close to pole")
	}
	a = compute(z, roots, &poles[1..]).map(|x| -x/(z_el.powf(2.0)));
	b = derive(z, roots, &poles[1..]).map(|x| x/z_el);
    } else {
	let z_el = z-roots[0];
	a = compute(z, &roots[1..], poles);
	b = derive(z, &roots[1..], poles).map(|x| x*z_el);
    };
    a.and_then(|x| b.map(|y| x+y))
}

#[pyfunction]
pub fn do_newton_iterate(maxiter: i64, z1: NewtComplex, roots: Vec<NewtComplex>,
			 poles: Vec<NewtComplex>) -> IterEnd {
    let mut i = 0_i64;
    let mut z = z1;
    let mut roots_copy = roots.to_vec();
    while i < maxiter {
	let derivative: NewtComplex;
	let value: NewtComplex;
	match derive(z, roots.as_slice(), poles.as_slice()) {
	    Ok(x) => derivative = x,
	    _ => break
	}

    	if derivative.norm_sqr() <= EPSILON_SQR.powf(2.0) {
	    break
	}

	match compute(z, roots.as_slice(), poles.as_slice()) {
	    Ok(x) => value = x,
	    _ => break
	}

	if value.norm_sqr() <= EPSILON_SQR {
	    roots_copy.sort_by(|k1, k2| (k1-z).norm_sqr().partial_cmp(&(k2-z).norm_sqr()).unwrap());
	    return IterEnd{
		z: Some(roots_copy[0]),
		niter: i,
	    }
	} else {
	    let step = -value/derivative;
	    z += step;
	    i += 1;
	}
    }
    IterEnd {
	z: None,
	niter: maxiter
    }
}










#[pymodule]
pub fn newton_iterate(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(do_newton_iterate))?;
    Ok(())
}
