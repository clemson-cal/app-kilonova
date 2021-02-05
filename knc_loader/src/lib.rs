use std::sync::Arc;
use pyo3::prelude::*;
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::PyIterProtocol;
use pyo3::PyMappingProtocol;
use pyo3::wrap_pyfunction;
use numpy::ToPyArray;
use pythonize::pythonize;
use kilonova::app;
use kilonova::io;
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
    products: Arc<products::Products>,
}

#[pyclass]
struct RadialProfileGetter {
    products: Arc<products::Products>,
}

#[pyclass]
struct RadialProfile {
    products: Arc<products::Products>,
    polar_index: usize,
}

#[pyclass]
struct BlockProducts {
    block_products: products::BlockProducts,
}




// ============================================================================
#[pymethods]
impl App {

    #[getter]
    fn version(&self, py: Python) -> PyResult<PyObject> {
        Ok(pythonize(py, &self.app.version)?)
    }

    #[getter]
    fn config(&self, py: Python) -> PyResult<PyObject> {
        Ok(pythonize(py, &self.app.config)?)
    }

    #[getter]
    fn tasks(&self, py: Python) -> PyResult<PyObject> {
        Ok(pythonize(py, &self.app.tasks)?)
    }

    fn make_products(&self) -> Products {
        Products{products: Arc::new(products::Products::from_app(&self.app))}
    }
}




// ============================================================================
#[pymethods]
impl Products {
    #[getter]
    fn time(&self) -> f64 {
        self.products.time
    }

    #[getter]
    fn config(&self, py: Python) -> PyResult<PyObject> {
        Ok(pythonize(py, &self.products.config)?)
    }

    #[getter]
    fn radial_profile(&self) -> RadialProfileGetter {
        RadialProfileGetter{products: self.products.clone()}
    }

    fn save(&self, filename: &str) -> PyResult<()> {
        match io::write_cbor(self.products.as_ref(), filename, false) {
            Ok(()) => Ok(()),
            Err(e) => Err(PyValueError::new_err(format!("{}", e))),
        }
    }
}




// ============================================================================
impl RadialProfile {

    fn concat_vertices(&self) -> ndarray::Array<f64, ndarray::Ix1> {
        let mut block_indexes: Vec<_> = self.products.blocks.keys().collect();
        block_indexes.sort();

        let arrays: Vec<_> = block_indexes
            .iter()
            .map(|i| self
                .products
                .blocks[i]
                .radial_vertices
                .view())
            .collect();

        ndarray::concatenate(ndarray::Axis(0), &arrays).unwrap()
    }

    fn concat_scalar(&self) -> ndarray::Array<f64, ndarray::Ix1> {
        let mut block_indexes: Vec<_> = self.products.blocks.keys().collect();
        block_indexes.sort();

        let arrays: Vec<_> = block_indexes
            .iter()
            .map(|i| self
                .products
                .blocks[i]
                .scalar
                .slice(ndarray::s![.., self.polar_index]))
            .collect();

        ndarray::concatenate(ndarray::Axis(0), &arrays).unwrap()
    }

    fn concat_map_primitive<F>(&self, f: F) -> ndarray::Array<f64, ndarray::Ix1>
    where
        F: Fn(&physics::AgnosticPrimitive) -> f64
    {
        let mut block_indexes: Vec<_> = self.products.blocks.keys().collect();
        block_indexes.sort();

        let arrays: Vec<_> = block_indexes
            .iter()
            .map(|i| self
                .products
                .blocks[i]
                .primitive
                .slice(ndarray::s![.., self.polar_index])
                .map(&f))
            .collect();
        let arrays: Vec<_> = arrays.iter().map(|a| a.view()).collect();

        ndarray::concatenate(ndarray::Axis(0), &arrays).unwrap()
    }
}

#[pymethods]
impl RadialProfile {

    #[getter]
    fn vertices(&self, py: Python) -> PyObject {
        self.concat_vertices().to_pyarray(py).to_object(py)
    }

    #[getter]
    fn scalar(&self, py: Python) -> PyObject {
        self.concat_scalar().to_pyarray(py).to_object(py)
    }

    #[getter]
    fn radial_four_velocity(&self, py: Python) -> PyObject {
        self.concat_map_primitive(|p| p.velocity_r).to_pyarray(py).to_object(py)
    }

    #[getter]
    fn polar_four_velocity(&self, py: Python) -> PyObject {
        self.concat_map_primitive(|p| p.velocity_q).to_pyarray(py).to_object(py)
    }

    #[getter]
    fn comoving_mass_density(&self, py: Python) -> PyObject {
        self.concat_map_primitive(|p| p.mass_density).to_pyarray(py).to_object(py)
    }

    #[getter]
    fn gas_pressure(&self, py: Python) -> PyObject {
        self.concat_map_primitive(|p| p.gas_pressure).to_pyarray(py).to_object(py)
    }
}




// ============================================================================
#[pyproto]
impl PyMappingProtocol for RadialProfileGetter {
    fn __getitem__(&self, polar_index: usize) -> RadialProfile {
        RadialProfile{products: self.products.clone(), polar_index}
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
            pyo3::Python::with_gil(|py| {
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
    fn map_primitive<F>(&self, f: F) -> ndarray::Array<f64, ndarray::Ix2>
    where
        F: Fn(&physics::AgnosticPrimitive) -> f64
    {
        self.block_products.primitive.map(f)
    }
}

#[pymethods]
impl BlockProducts {

    #[getter]
    fn radial_vertices(&self, py: Python) -> PyObject {
        self.block_products.radial_vertices.to_pyarray(py).to_object(py)
    }

    #[getter]
    fn polar_vertices(&self, py: Python) -> PyObject {
        self.block_products.polar_vertices.to_pyarray(py).to_object(py)
    }

    #[getter]
    fn scalar(&self, py: Python) -> PyObject {
        self.block_products.scalar.to_pyarray(py).to_object(py)
    }

    #[getter]
    fn radial_four_velocity(&self, py: Python) -> PyObject {
        self.map_primitive(|p| p.velocity_r).to_pyarray(py).to_object(py)
    }

    #[getter]
    fn polar_four_velocity(&self, py: Python) -> PyObject {
        self.map_primitive(|p| p.velocity_q).to_pyarray(py).to_object(py)
    }

    #[getter]
    fn comoving_mass_density(&self, py: Python) -> PyObject {
        self.map_primitive(|p| p.mass_density).to_pyarray(py).to_object(py)
    }

    #[getter]
    fn gas_pressure(&self, py: Python) -> PyObject {
        self.map_primitive(|p| p.gas_pressure).to_pyarray(py).to_object(py)
    }
}




// ============================================================================
#[pyfunction]
fn app(filename: &str) -> PyResult<App> {
    match app::App::from_preset_or_file(filename) {
        Ok(app) => Ok(App{app}),
        Err(e)  => Err(PyValueError::new_err(format!("{}", e))),
    }
}

#[pyfunction]
fn products(filename: &str) -> PyResult<Products> {
    match io::read_cbor(filename, false) {
        Ok(products) => Ok(Products{products: Arc::new(products)}),
        Err(e)       => Err(PyValueError::new_err(format!("{}", e))),
    }
}




// ============================================================================
#[pymodule]
fn knc_loader(_: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(app, m)?)?;
    m.add_function(wrap_pyfunction!(products, m)?)?;
    Ok(())
}
