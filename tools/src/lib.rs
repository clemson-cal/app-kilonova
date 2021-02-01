use pyo3::prelude::*;
use pyo3::exceptions::PyKeyError;
use pyo3::PyIterProtocol;
use pyo3::PyMappingProtocol;
use pyo3::wrap_pyfunction;
use numpy::ToPyArray;
use pythonize::pythonize;
use kilonova::app;
use kilonova::mesh;
use kilonova::physics;
use kilonova::products;




// ============================================================================
#[pyclass]
struct App {
	app: app::App
}

#[pyclass]
struct Products {
	products: products::Products,
}

#[pyclass]
struct BlockProducts {
	block_products: products::BlockProducts,
}




// ============================================================================
#[pymethods]
impl App {

	#[getter]
	fn version(&self) -> PyResult<PyObject> {
		Python::with_gil(|py| {
			Ok(pythonize(py, &self.app.version)?)
		})
	}

	#[getter]
	fn config(&self) -> PyResult<PyObject> {
		Python::with_gil(|py| {
			Ok(pythonize(py, &self.app.config)?)
		})
	}

	#[getter]
	fn tasks(&self) -> PyResult<PyObject> {
		Python::with_gil(|py| {
			Ok(pythonize(py, &self.app.tasks)?)
		})
	}

	fn make_products(&self) -> Products {
		Products{products: products::Products::from_app(&self.app)}
	}
}




// ============================================================================
#[pymethods]
impl Products {
	#[getter]
	fn time(&self) -> f64 {
		self.products.time
	}
}




// ============================================================================
#[pyproto]
impl PyMappingProtocol for Products {

	fn __len__(&self) -> usize {
		self.products.blocks.len()
	}

	fn __getitem__(&self, key: mesh::BlockIndex) -> PyResult<BlockProducts> {
		if let Some(b) = self.products.blocks.get(&key) {
			Ok(BlockProducts{block_products: b.clone()})
		} else {
			Python::with_gil(|py| {
				Err(PyErr::from_instance(PyKeyError::new_err("invalid block index").instance(py)))
			})
		}
	}
}

#[pyproto]
impl PyIterProtocol for Products {
    fn __iter__(slf: PyRef<Self>) -> PyResult<Py<ProductsIter>> {
    	let keys: Vec<_> = slf.products.blocks.keys().cloned().collect();
        let iter = ProductsIter {
            inner: keys.into_iter()
        };
        Py::new(slf.py(), iter)
    }
}




// ============================================================================
#[pyclass]
struct ProductsIter {
    inner: std::vec::IntoIter<mesh::BlockIndex>,
}

#[pyproto]
impl PyIterProtocol for ProductsIter {

    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<Self>) -> Option<mesh::BlockIndex> {
        slf.inner.next()
    }
}




// ============================================================================
impl BlockProducts {
	fn map_primitive<F>(&self, f: F) -> PyObject
	where
		F: Fn(&physics::AgnosticPrimitive) -> f64
	{
		pyo3::Python::with_gil(|py| {
			self.block_products.primitive.map(f).to_pyarray(py).to_object(py)
		})
	}
}

#[pymethods]
impl BlockProducts {

	#[getter]
	fn radial_vertices(&self) -> PyObject {
		pyo3::Python::with_gil(|py| {
			self.block_products.radial_vertices.to_pyarray(py).to_object(py)
		})
	}

	#[getter]
	fn polar_vertices(&self) -> PyObject {
		pyo3::Python::with_gil(|py| {
			self.block_products.polar_vertices.to_pyarray(py).to_object(py)
		})
	}

	#[getter]
	fn scalar(&self) -> PyObject {
		pyo3::Python::with_gil(|py| {
			self.block_products.scalar.to_pyarray(py).to_object(py)
		})
	}

	#[getter]
	fn radial_four_velocity(&self) -> PyObject {
		self.map_primitive(|p| p.velocity_r)
	}

	#[getter]
	fn polar_four_velocity(&self) -> PyObject {
		self.map_primitive(|p| p.velocity_q)
	}

	#[getter]
	fn comoving_mass_density(&self) -> PyObject {
		self.map_primitive(|p| p.mass_density)
	}

	#[getter]
	fn gas_pressure(&self) -> PyObject {
		self.map_primitive(|p| p.gas_pressure)
	}
}




// ============================================================================
#[pyfunction]
fn app(filename: &str) -> App {
	let app = app::App::from_preset_or_file(filename).unwrap();
	App{app}
}




// ============================================================================
#[pymodule]
fn knc_tools(_: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(app, m)?)?;
    Ok(())
}
